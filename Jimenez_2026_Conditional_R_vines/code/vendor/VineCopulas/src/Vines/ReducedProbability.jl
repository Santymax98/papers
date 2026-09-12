# Experimental reduced probability integrals for full simplified continuous
# vines. Everything in this file is private: no public cdf method dispatches
# here while the numerical research prototype is being evaluated.

struct _TerminalPairInfo{C<:PairCopula}
    a::Int
    b::Int
    conditioning::Vector{Int}
    copula::C
    tree::Int
    index::Int
end

struct _ConditioningSubvineDescriptor
    nodes::Vector{Int}
    edges::Vector{VineEdge}
end

struct _ReducedStateOp
    left::Int
    right::Int
    out_left::Int
    out_right::Int
    copula_index::Int
    tree::Int
    edge::Int
end

struct _ReducedProbabilityPlan{T,C,PC}
    p::Int
    terminal::T
    terminal_copula::C
    copulas::PC
    subvine::_ConditioningSubvineDescriptor
    rawslots::Vector{Int}
    states::Dict{Any,Int}
    slot_variable::Vector{Int}
    slot_conditioning::Vector{Vector{Int}}
    invalidate_slots::Vector{Vector{Int}}
    producer::Vector{Int}
    producer_side::Vector{UInt8}
    ops::Vector{_ReducedStateOp}
    default_order::Vector{Int}
    terminal_a_slot::Int
    terminal_b_slot::Int
end

mutable struct _ReducedEvalCounts
    hfunc::Int
    hinv::Int
    pair_cdf::Int
end

_ReducedEvalCounts() = _ReducedEvalCounts(0, 0, 0)

# The original vine stores pair copulas in a concrete tuple. The execution
# plan records integer copula indices in state ops and stores homogeneous vines
# in a concrete vector for O(1) indexed access; heterogeneous vines retain the
# concrete tuple used by the generated switch below.
@generated function _reduced_copula_storage(copulas::PC) where {PC<:Tuple}
    types = fieldtypes(PC)
    return all(==(first(types)), types) ? :(collect(copulas)) : :copulas
end

@generated function _reduced_apply_pair(
    ::Val{F}, copulas::PC, index::Int, x::Float64, y::Float64,
)::Float64 where {F,PC}
    apply = if F === :hfunc1
        :hfunc1
    elseif F === :hfunc2
        :hfunc2
    elseif F === :hinv1
        :hinv1
    elseif F === :hinv2
        :hinv2
    else
        error("unknown reduced pair operation $F")
    end
    if PC <: AbstractVector
        return quote
            @inbounds return $apply(copulas[index], x, y)::Float64
        end
    end
    PC <: Tuple || error("unsupported compiled pair-copula storage $PC")
    cases = Expr(:block)
    for i in 1:fieldcount(PC)
        push!(cases.args, quote
            if index == $i
                return $apply(getfield(copulas, $i), x, y)::Float64
            end
        end)
    end
    return quote
        @inbounds begin
            $cases
        end
        throw(BoundsError(copulas, index))
    end
end

mutable struct _ReducedProbabilityWorkspace
    values::Vector{Float64}
    valid::BitVector
    cdfbuf::Vector{Float64}
    counts::_ReducedEvalCounts
end

struct _ReducedRQMCResult
    estimate::Float64
    stderr::Float64
    replicates::Vector{Float64}
    N::Int
    R::Int
    dimension::Int
    order::Vector{Int}
    hfunc_calls::Int
    hinv_calls::Int
    pair_cdf_calls::Int
end

@inline _reduced_state_key(v::Int, D) = (v, Tuple(sort!(collect(Int, D))))

function _terminal_pair_info(vc::AbstractVineCopula{p}) where {p}
    p >= 3 || throw(ArgumentError("reduced probability integration requires dimension at least 3"))
    truncation(vc) == p - 1 || throw(ArgumentError(
        "reduced probability integration currently requires a full, non-truncated vine",
    ))
    terminal_edges = filter(e -> e.tree == p - 1, vine_edges(vc))
    length(terminal_edges) == 1 || throw(ArgumentError(
        "a full vine must contain exactly one edge in its last tree; found $(length(terminal_edges))",
    ))
    edge = only(terminal_edges)
    a, b = edge.conditioned
    D = collect(Int, edge.conditioning)
    length(D) == p - 2 || throw(ArgumentError("terminal conditioning set has wrong dimension"))
    length(unique(D)) == length(D) || throw(ArgumentError("terminal conditioning set has duplicates"))
    return _TerminalPairInfo(a, b, D, edge.copula, edge.tree, edge.index)
end

function _conditioning_subvine_descriptor(vc::AbstractVineCopula)
    terminal = _terminal_pair_info(vc)
    Dset = Set(terminal.conditioning)
    subedges = VineEdge[]
    for edge in vine_edges(vc)
        support = Set((edge.conditioned..., edge.conditioning...))
        support <= Dset && push!(subedges, edge)
    end
    s = length(terminal.conditioning)
    expected = s * (s - 1) ÷ 2
    length(subedges) == expected || throw(ArgumentError(
        "terminal conditioning set does not induce a complete sub-vine: expected $expected edges, found $(length(subedges))",
    ))
    sort!(subedges; by=e -> (e.tree, e.index))
    return _ConditioningSubvineDescriptor(sort(copy(terminal.conditioning)), subedges)
end

function _compile_reduced_probability(vc::AbstractVineCopula{p}) where {p}
    terminal = _terminal_pair_info(vc)
    original_copulas = Tuple(Iterators.flatten(vc.edges))
    terminal_copula = last(original_copulas)
    terminal.copula === terminal_copula || throw(ArgumentError(
        "terminal edge does not match the final pair copula stored by the vine",
    ))
    _pair_cdf_is_safe(terminal_copula) || throw(ArgumentError(
        "terminal pair CDF for $(typeof(terminal.copula)) is not classified as deterministic; " *
        "outer RQMC refuses an unknown or potentially randomized inner CDF",
    ))
    subvine = _conditioning_subvine_descriptor(vc)
    copulas = _reduced_copula_storage(original_copulas)

    states = Dict{Any,Int}()
    rawslots = zeros(Int, p)
    slot_variable = Int[]
    slot_conditioning = Vector{Int}[]
    producer = Int[]
    producer_side = UInt8[]
    for v in 1:p
        slot = length(slot_variable) + 1
        states[_reduced_state_key(v, Int[])] = slot
        rawslots[v] = slot
        push!(slot_variable, v)
        push!(slot_conditioning, Int[])
        push!(producer, 0)
        push!(producer_side, 0x00)
    end

    ops = _ReducedStateOp[]
    all_edges = sort(vine_edges(vc); by=e -> (e.tree, e.index))
    length(all_edges) == length(original_copulas) || throw(ArgumentError(
        "vine edge metadata and stored pair-copula sequence have different lengths",
    ))
    for (copula_index, edge) in enumerate(all_edges)
        a, b = edge.conditioned
        D = collect(Int, edge.conditioning)
        ka, kb = _reduced_state_key(a, D), _reduced_state_key(b, D)
        haskey(states, ka) || throw(ArgumentError("missing conditional state $ka"))
        haskey(states, kb) || throw(ArgumentError("missing conditional state $kb"))
        oa = _reduced_state_key(a, vcat(D, b))
        ob = _reduced_state_key(b, vcat(D, a))
        (haskey(states, oa) || haskey(states, ob)) && throw(ArgumentError(
            "conditional state generated more than once by vine metadata",
        ))
        sa = length(slot_variable) + 1
        sb = sa + 1
        states[oa], states[ob] = sa, sb
        opidx = length(ops) + 1
        push!(ops, _ReducedStateOp(
            states[ka], states[kb], sa, sb, copula_index, edge.tree, edge.index,
        ))
        for (v, cond, side) in ((a, sort!(vcat(copy(D), b)), 0x01),
                                (b, sort!(vcat(copy(D), a)), 0x02))
            push!(slot_variable, v)
            push!(slot_conditioning, cond)
            push!(producer, opidx)
            push!(producer_side, side)
        end
    end

    # The conditioning-set order stored by the terminal edge records the
    # natural construction order for D/C/R metadata. Validate it through the
    # same general scheduler rather than trusting a family convention.
    provisional = _ReducedProbabilityPlan(
        p, terminal, terminal_copula, copulas, subvine, rawslots, states, slot_variable,
        slot_conditioning, Vector{Vector{Int}}(), producer, producer_side,
        ops, Int[], 0, 0,
    )
    default_order = _reduced_sampling_order(provisional, zeros(p), :default)
    terminal_a_slot = states[_reduced_state_key(terminal.a, terminal.conditioning)]
    terminal_b_slot = states[_reduced_state_key(terminal.b, terminal.conditioning)]
    invalidate_slots = [Int[] for _ in 1:p]
    for variable in 1:p, slot in eachindex(slot_variable)
        if slot_variable[slot] == variable || variable in slot_conditioning[slot]
            push!(invalidate_slots[variable], slot)
        end
    end
    return _ReducedProbabilityPlan(
        p, terminal, terminal_copula, copulas, subvine, rawslots, states, slot_variable,
        slot_conditioning, invalidate_slots, producer, producer_side, ops, default_order,
        terminal_a_slot, terminal_b_slot,
    )
end

function _reduced_order_search!(
    plan::_ReducedProbabilityPlan,
    prefix::Vector{Int},
    remaining::Vector{Int},
    preference::Vector{Int},
)
    isempty(remaining) && return copy(prefix)
    Dkey = sort(copy(prefix))
    for candidate in preference
        candidate in remaining || continue
        haskey(plan.states, _reduced_state_key(candidate, Dkey)) || continue
        push!(prefix, candidate)
        nextremaining = filter(!=(candidate), remaining)
        result = _reduced_order_search!(plan, prefix, nextremaining, preference)
        result === nothing || return result
        pop!(prefix)
    end
    return nothing
end

function _reduced_sampling_order(
    plan::_ReducedProbabilityPlan,
    widths::AbstractVector{<:Real},
    strategy::Symbol,
)
    strategy in (:default, :smallest_width_first) || throw(ArgumentError(
        "order_strategy must be :default or :smallest_width_first",
    ))
    D = copy(plan.terminal.conditioning)
    preference = if strategy === :default
        D
    else
        sort(copy(D); by=v -> (Float64(widths[v]), findfirst(==(v), D)))
    end
    result = _reduced_order_search!(plan, Int[], copy(D), preference)
    result === nothing && throw(ArgumentError(
        "could not construct an admissible sampling order for terminal conditioning set $D",
    ))
    return result
end

function _reduced_workspace(plan::_ReducedProbabilityPlan)
    n = length(plan.slot_variable)
    return _ReducedProbabilityWorkspace(zeros(n), falses(n), zeros(2), _ReducedEvalCounts())
end

@inline function _invalidate_variable!(
    ws::_ReducedProbabilityWorkspace,
    plan::_ReducedProbabilityPlan,
    variable::Int,
)
    @inbounds for slot in plan.invalidate_slots[variable]
        ws.valid[slot] = false
    end
    return nothing
end

@inline function _set_raw_value!(
    ws::_ReducedProbabilityWorkspace,
    plan::_ReducedProbabilityPlan,
    variable::Int,
    value::Real,
)
    _invalidate_variable!(ws, plan, variable)
    slot = plan.rawslots[variable]
    ws.values[slot] = clamp(Float64(value), 0.0, 1.0)
    ws.valid[slot] = true
    return ws.values[slot]
end

function _eval_reduced_state!(
    ws::_ReducedProbabilityWorkspace,
    plan::_ReducedProbabilityPlan,
    slot::Int,
)::Float64
    ws.valid[slot] && return ws.values[slot]
    opidx = plan.producer[slot]
    opidx > 0 || throw(ArgumentError("raw conditional state has no assigned value"))
    op = plan.ops[opidx]
    left = _eval_reduced_state!(ws, plan, op.left)
    right = _eval_reduced_state!(ws, plan, op.right)
    value = if plan.producer_side[slot] == 0x01
        _reduced_apply_pair(Val(:hfunc1), plan.copulas, op.copula_index, left, right)
    else
        _reduced_apply_pair(Val(:hfunc2), plan.copulas, op.copula_index, left, right)
    end
    ws.counts.hfunc += 1
    value = clamp(Float64(value), 0.0, 1.0)
    ws.values[slot] = value
    ws.valid[slot] = true
    return value
end

function _conditional_endpoint!(
    ws::_ReducedProbabilityWorkspace,
    plan::_ReducedProbabilityPlan,
    variable::Int,
    conditioning,
    endpoint::Real,
)::Float64
    x = clamp(Float64(endpoint), 0.0, 1.0)
    iszero(x) && return 0.0
    isone(x) && return 1.0
    _set_raw_value!(ws, plan, variable, x)
    key = _reduced_state_key(variable, conditioning)
    haskey(plan.states, key) || throw(ArgumentError("conditional state $key is unavailable"))
    return _eval_reduced_state!(ws, plan, plan.states[key])
end

@inline function _conditional_endpoint_slot!(
    ws::_ReducedProbabilityWorkspace,
    plan::_ReducedProbabilityPlan,
    variable::Int,
    slot::Int,
    endpoint::Real,
)::Float64
    x = clamp(Float64(endpoint), 0.0, 1.0)
    iszero(x) && return 0.0
    isone(x) && return 1.0
    _set_raw_value!(ws, plan, variable, x)
    return _eval_reduced_state!(ws, plan, slot)
end

function _reduced_sampling_slots(
    plan::_ReducedProbabilityPlan,
    sampling_order::AbstractVector{<:Integer},
)
    slots = Vector{Int}(undef, length(sampling_order))
    prefix = Int[]
    for (j, variable0) in enumerate(sampling_order)
        variable = Int(variable0)
        slots[j] = plan.states[_reduced_state_key(variable, prefix)]
        push!(prefix, variable)
    end
    return slots
end

function _invert_reduced_state!(
    ws::_ReducedProbabilityWorkspace,
    plan::_ReducedProbabilityPlan,
    slot::Int,
    probability::Real,
)::Float64
    opidx = plan.producer[slot]
    opidx == 0 && return clamp(Float64(probability), 0.0, 1.0)
    op = plan.ops[opidx]
    q = clamp(Float64(probability), 0.0, 1.0)
    if plan.producer_side[slot] == 0x01
        other = _eval_reduced_state!(ws, plan, op.right)
        parent, value = op.left, _reduced_apply_pair(
            Val(:hinv1), plan.copulas, op.copula_index, q, other,
        )
    else
        other = _eval_reduced_state!(ws, plan, op.left)
        parent, value = op.right, _reduced_apply_pair(
            Val(:hinv2), plan.copulas, op.copula_index, q, other,
        )
    end
    ws.counts.hinv += 1
    return _invert_reduced_state!(ws, plan, parent, value)
end

function _pair_rectprob_counted!(
    ws::_ReducedProbabilityWorkspace,
    C::PairCopula,
    lo1::Float64,
    hi1::Float64,
    lo2::Float64,
    hi2::Float64,
)::Float64
    (hi1 <= lo1 || hi2 <= lo2) && return 0.0
    p = _pair_cdf(C, hi1, hi2, ws.cdfbuf)
    ws.counts.pair_cdf += 1
    if !iszero(lo1)
        p -= _pair_cdf(C, lo1, hi2, ws.cdfbuf)
        ws.counts.pair_cdf += 1
    end
    if !iszero(lo2)
        p -= _pair_cdf(C, hi1, lo2, ws.cdfbuf)
        ws.counts.pair_cdf += 1
    end
    if !iszero(lo1) && !iszero(lo2)
        p += _pair_cdf(C, lo1, lo2, ws.cdfbuf)
        ws.counts.pair_cdf += 1
    end
    return clamp(Float64(p), 0.0, 1.0)
end

function _truncated_inverse_rosenblatt!(
    values::AbstractVector{Float64},
    ws::_ReducedProbabilityWorkspace,
    plan::_ReducedProbabilityPlan,
    lower::AbstractVector{<:Real},
    upper::AbstractVector{<:Real},
    z::AbstractVector{<:Real},
    sampling_order::AbstractVector{<:Integer},
    sampling_slots::AbstractVector{<:Integer},
)::Tuple{Float64,Float64}
    fill!(ws.valid, false)
    fill!(values, NaN)
    weight = 1.0
    q1 = 1.0
    @inbounds for (j, variable0) in enumerate(sampling_order)
        variable = Int(variable0)
        alpha = if j == 1
            Float64(lower[variable])
        else
            _conditional_endpoint_slot!(
                ws, plan, variable, sampling_slots[j], lower[variable],
            )
        end
        beta = if j == 1
            Float64(upper[variable])
        else
            _conditional_endpoint_slot!(
                ws, plan, variable, sampling_slots[j], upper[variable],
            )
        end
        delta = beta - alpha
        delta >= -256eps(Float64) || throw(DomainError(
            delta, "conditional interval has negative probability; check h-function orientation",
        ))
        q = clamp(delta, 0.0, 1.0)
        j == 1 && (q1 = q)
        weight *= q
        iszero(weight) && return 0.0, q1
        r = clamp(alpha + clamp(Float64(z[j]), 0.0, 1.0) * q, alpha, beta)
        state = sampling_slots[j]
        value = j == 1 ? r : _invert_reduced_state!(ws, plan, state, r)
        values[variable] = _set_raw_value!(ws, plan, variable, value)
    end
    return weight, q1
end

function _reduced_probability_integrand!(
    ws::_ReducedProbabilityWorkspace,
    plan::_ReducedProbabilityPlan,
    lower::AbstractVector{<:Real},
    upper::AbstractVector{<:Real},
    z::AbstractVector{<:Real},
    sampling_order::AbstractVector{<:Integer},
    sampling_slots::AbstractVector{<:Integer},
    conditioning_values::AbstractVector{Float64},
)::Tuple{Float64,Float64}
    weight, q1 = _truncated_inverse_rosenblatt!(
        conditioning_values, ws, plan, lower, upper, z, sampling_order,
        sampling_slots,
    )
    iszero(weight) && return 0.0, q1
    terminal = plan.terminal
    D = terminal.conditioning
    aminus = _conditional_endpoint_slot!(
        ws, plan, terminal.a, plan.terminal_a_slot, lower[terminal.a],
    )
    aplus = _conditional_endpoint_slot!(
        ws, plan, terminal.a, plan.terminal_a_slot, upper[terminal.a],
    )
    bminus = _conditional_endpoint_slot!(
        ws, plan, terminal.b, plan.terminal_b_slot, lower[terminal.b],
    )
    bplus = _conditional_endpoint_slot!(
        ws, plan, terminal.b, plan.terminal_b_slot, upper[terminal.b],
    )
    H = _pair_rectprob_counted!(
        ws, plan.terminal_copula, aminus, aplus, bminus, bplus,
    )
    G = weight * H
    tolerance = 2048eps(Float64) * max(1.0, q1)
    (-tolerance <= G <= q1 + tolerance) || throw(DomainError(
        G, "reduced integrand violates 0 <= G <= q1; check conditional orientation",
    ))
    return clamp(G, 0.0, q1), q1
end

function _check_reduced_bounds(
    p::Int,
    lower::AbstractVector{<:Real},
    upper::AbstractVector{<:Real},
)
    _check_vector_dim(p, lower)
    _check_vector_dim(p, upper)
    lo, hi = Vector{Float64}(undef, p), Vector{Float64}(undef, p)
    @inbounds for i in 1:p
        l, u = Float64(lower[i]), Float64(upper[i])
        (isfinite(l) && isfinite(u)) || throw(DomainError((l, u), "bounds must be finite"))
        (0.0 <= l <= u <= 1.0) || throw(DomainError((l, u), "bounds must satisfy 0 <= lower <= upper <= 1"))
        lo[i], hi[i] = l, u
    end
    return lo, hi
end

function _reduced_owen_points(dimension::Int, N::Int, seed::UInt64; randomized::Bool=true)
    N >= 1 || throw(ArgumentError("N must be positive"))
    ispow2(N) || throw(ArgumentError("experimental reduced RQMC requires N=2^m"))
    X = _qmc_points(dimension, N; randomized=false)
    randomized || return X
    Y = similar(X)
    @inbounds for axis in 1:dimension
        axis_seed = seed + UInt64(axis) * 0x9e3779b97f4a7c15
        rng = Random.MersenneTwister(axis_seed)
        row = Matrix(reshape(copy(@view(X[axis, :])), 1, N))
        scrambled = QuasiMonteCarlo.randomize(
            row, QuasiMonteCarlo.OwenScramble(base=2, rng=rng),
        )
        @views Y[axis, :] .= vec(scrambled)
    end
    return Y
end

function _rectprob_reduced_rqmc(
    vc::AbstractVineCopula{p},
    lower::AbstractVector{<:Real},
    upper::AbstractVector{<:Real};
    N::Integer=1 << 12,
    R::Integer=8,
    seed::Integer=0x5245514d,
    order_strategy::Symbol=:default,
    randomized::Bool=true,
) where {p}
    Ni, Ri = Int(N), Int(R)
    Ri >= 1 || throw(ArgumentError("R must be positive"))
    lo, hi = _check_reduced_bounds(p, lower, upper)
    plan = _compile_reduced_probability(vc)
    sampling_order = _reduced_sampling_order(plan, hi .- lo, order_strategy)
    sampling_slots = _reduced_sampling_slots(plan, sampling_order)
    dimension = p - 2
    replicates = zeros(Float64, Ri)
    if any(i -> lo[i] == hi[i], 1:p)
        return _ReducedRQMCResult(0.0, 0.0, replicates, Ni, Ri, dimension,
                                  sampling_order, 0, 0, 0)
    end
    if all(iszero, lo) && all(isone, hi)
        fill!(replicates, 1.0)
        return _ReducedRQMCResult(1.0, 0.0, replicates, Ni, Ri, dimension,
                                  sampling_order, 0, 0, 0)
    end

    ws = _reduced_workspace(plan)
    conditioning_values = fill(NaN, p)
    for rep in 1:Ri
        Z = _reduced_owen_points(
            dimension, Ni, UInt64(seed) + UInt64(rep - 1); randomized=randomized,
        )
        total = 0.0
        @inbounds for col in axes(Z, 2)
            G, _ = _reduced_probability_integrand!(
                ws, plan, lo, hi, view(Z, :, col), sampling_order,
                sampling_slots, conditioning_values,
            )
            total += G
        end
        replicates[rep] = total / Ni
    end
    estimate = sum(replicates) / Ri
    stderr = Ri > 1 ? sqrt(sum(abs2, replicates .- estimate) / (Ri - 1)) / sqrt(Ri) : NaN
    return _ReducedRQMCResult(
        estimate, stderr, replicates, Ni, Ri, dimension, sampling_order,
        ws.counts.hfunc, ws.counts.hinv, ws.counts.pair_cdf,
    )
end

function _cdf_reduced_rqmc(
    vc::AbstractVineCopula{p},
    upper::AbstractVector{<:Real};
    kwargs...,
) where {p}
    return _rectprob_reduced_rqmc(vc, zeros(p), upper; kwargs...)
end

function _reduced_probability_integrand(
    vc::AbstractVineCopula{p},
    lower::AbstractVector{<:Real},
    upper::AbstractVector{<:Real},
    z::AbstractVector{<:Real};
    order_strategy::Symbol=:default,
) where {p}
    lo, hi = _check_reduced_bounds(p, lower, upper)
    plan = _compile_reduced_probability(vc)
    sampling_order = _reduced_sampling_order(plan, hi .- lo, order_strategy)
    sampling_slots = _reduced_sampling_slots(plan, sampling_order)
    length(z) == p - 2 || throw(DimensionMismatch("z must have length $(p-2)"))
    ws = _reduced_workspace(plan)
    values = fill(NaN, p)
    G, q1 = _reduced_probability_integrand!(
        ws, plan, lo, hi, z, sampling_order, sampling_slots, values,
    )
    return (;
        value=G,
        q1,
        order=sampling_order,
        conditioning_values=values,
        hfunc_calls=ws.counts.hfunc,
        hinv_calls=ws.counts.hinv,
        pair_cdf_calls=ws.counts.pair_cdf,
    )
end
