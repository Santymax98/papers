# Private research extensions for sequential probability regions and ordinal
# patterns. These routines deliberately do not extend the public cdf API.

function _reduced_admissible_orders(plan::_ReducedProbabilityPlan)
    results = Vector{Vector{Int}}()
    prefix = Int[]
    remaining = copy(plan.terminal.conditioning)
    function visit!()
        if isempty(remaining)
            push!(results, copy(prefix))
            return
        end
        prefix_key = sort(copy(prefix))
        for variable in copy(remaining)
            haskey(plan.states, _reduced_state_key(variable, prefix_key)) || continue
            index = findfirst(==(variable), remaining)
            deleteat!(remaining, index)
            push!(prefix, variable)
            visit!()
            pop!(prefix)
            insert!(remaining, index, variable)
        end
    end
    visit!()
    sort!(results)
    return results
end

function _validate_reduced_order(
    plan::_ReducedProbabilityPlan,
    sampling_order::AbstractVector{<:Integer},
)
    order = collect(Int, sampling_order)
    sort(order) == sort(plan.terminal.conditioning) || throw(ArgumentError(
        "sampling order must be a permutation of the terminal conditioning set",
    ))
    prefix = Int[]
    for variable in order
        haskey(plan.states, _reduced_state_key(variable, prefix)) || throw(ArgumentError(
            "sampling order $(join(order, '-')) is not reachable through the vine state recursion",
        ))
        push!(prefix, variable)
    end
    return order
end

"""
    _truncated_inverse_rosenblatt_region!(..., bounds!)

Private generic sequential-region transport. `bounds!(variable, stage,
values, sampled)` returns the raw-scale lower and upper bounds for the next
variable. The routine performs only the same conditional-CDF and h-inverse
operations used by the rectangular engine.
"""
function _truncated_inverse_rosenblatt_region!(
    values::AbstractVector{Float64},
    sampled::AbstractVector{Bool},
    ws::_ReducedProbabilityWorkspace,
    plan::_ReducedProbabilityPlan,
    z::AbstractVector{<:Real},
    sampling_order::AbstractVector{<:Integer},
    sampling_slots::AbstractVector{<:Integer},
    bounds!,
)::Tuple{Float64,Float64}
    fill!(ws.valid, false)
    fill!(values, NaN)
    fill!(sampled, false)
    weight = 1.0
    q1 = 1.0
    @inbounds for (stage, variable0) in enumerate(sampling_order)
        variable = Int(variable0)
        lower, upper = bounds!(variable, stage, values, sampled)
        0.0 <= lower <= upper <= 1.0 || throw(DomainError(
            (lower, upper), "sequential region bounds must lie in [0,1] and be ordered",
        ))
        alpha = stage == 1 ? Float64(lower) : _conditional_endpoint_slot!(
            ws, plan, variable, sampling_slots[stage], lower,
        )
        beta = stage == 1 ? Float64(upper) : _conditional_endpoint_slot!(
            ws, plan, variable, sampling_slots[stage], upper,
        )
        delta = beta - alpha
        delta >= -256eps(Float64) || throw(DomainError(
            delta, "conditional sequential interval has negative probability",
        ))
        q = clamp(delta, 0.0, 1.0)
        stage == 1 && (q1 = q)
        weight *= q
        iszero(weight) && return 0.0, q1
        r = clamp(alpha + clamp(Float64(z[stage]), 0.0, 1.0) * q, alpha, beta)
        state = sampling_slots[stage]
        value = stage == 1 ? r : _invert_reduced_state!(ws, plan, state, r)
        values[variable] = _set_raw_value!(ws, plan, variable, value)
        sampled[variable] = true
    end
    return weight, q1
end

function _ordinal_ranks(pattern::AbstractVector{<:Integer}, p::Int)
    length(pattern) == p || throw(DimensionMismatch("ordinal pattern must have length $p"))
    sort(collect(Int, pattern)) == collect(1:p) || throw(ArgumentError(
        "ordinal pattern must be a permutation of 1:$p",
    ))
    ranks = zeros(Int, p)
    @inbounds for (rank, variable0) in enumerate(pattern)
        ranks[Int(variable0)] = rank
    end
    return ranks
end

"""
    _normalized_permutation_entropy_with_gradient(raw_probabilities)

Private numerical helper for the JMVA application.  It returns normalized
permutation entropy and its gradient with respect to the unnormalized raw
probabilities.  The gradient accounts for the dependence of the normalizing
mass on every probability.
"""
function _normalized_permutation_entropy_with_gradient(
    raw_probabilities::AbstractVector{<:Real},
)
    K = length(raw_probabilities)
    K >= 2 || throw(ArgumentError("at least two probabilities are required"))
    mass = sum(raw_probabilities)
    isfinite(mass) && mass > 0 || throw(DomainError(
        mass, "raw probability mass must be finite and positive",
    ))
    probabilities = Float64.(raw_probabilities) ./ mass
    all(isfinite, probabilities) || throw(DomainError(
        probabilities, "probabilities must be finite",
    ))
    all(>(0.0), probabilities) || throw(DomainError(
        probabilities,
        "the ordinary delta-method gradient is undefined at zero probability",
    ))
    logK = log(K)
    entropy = -sum(p * log(p) for p in probabilities) / logK
    probability_gradient = [-(log(p) + 1) / logK for p in probabilities]
    centered = sum(probabilities .* probability_gradient)
    raw_gradient = (probability_gradient .- centered) ./ mass
    return (; entropy, gradient=raw_gradient, probabilities, mass)
end

@inline function _ordinal_gap(
    variable::Int,
    ranks::AbstractVector{<:Integer},
    values::AbstractVector{Float64},
    sampled::AbstractVector{Bool},
)
    target_rank = ranks[variable]
    lower_rank = 0
    upper_rank = length(ranks) + 1
    lower = 0.0
    upper = 1.0
    @inbounds for other in eachindex(ranks)
        sampled[other] || continue
        rank = ranks[other]
        if lower_rank < rank < target_rank
            lower_rank = rank
            lower = values[other]
        elseif target_rank < rank < upper_rank
            upper_rank = rank
            upper = values[other]
        end
    end
    return lower, upper
end

function _ordinal_probability_integrand!(
    ws::_ReducedProbabilityWorkspace,
    plan::_ReducedProbabilityPlan,
    pattern::AbstractVector{<:Integer},
    ranks::AbstractVector{<:Integer},
    z::AbstractVector{<:Real},
    sampling_order::AbstractVector{<:Integer},
    sampling_slots::AbstractVector{<:Integer},
    conditioning_values::AbstractVector{Float64},
    sampled::AbstractVector{Bool},
)::Tuple{Float64,Float64}
    terminal = plan.terminal
    adjacent = abs(ranks[terminal.a] - ranks[terminal.b]) == 1
    expected_dimension = plan.p - 2 + Int(adjacent)
    length(z) == expected_dimension || throw(DimensionMismatch(
        "ordinal integrand requires $expected_dimension coordinates",
    ))
    bounds! = (variable, stage, values, active) ->
        _ordinal_gap(variable, ranks, values, active)
    weight, q1 = _truncated_inverse_rosenblatt_region!(
        conditioning_values, sampled, ws, plan, view(z, 1:(plan.p - 2)),
        sampling_order, sampling_slots, bounds!,
    )
    iszero(weight) && return 0.0, q1

    if !adjacent
        aminus_raw, aplus_raw = _ordinal_gap(
            terminal.a, ranks, conditioning_values, sampled,
        )
        bminus_raw, bplus_raw = _ordinal_gap(
            terminal.b, ranks, conditioning_values, sampled,
        )
        aminus = _conditional_endpoint_slot!(
            ws, plan, terminal.a, plan.terminal_a_slot, aminus_raw,
        )
        aplus = _conditional_endpoint_slot!(
            ws, plan, terminal.a, plan.terminal_a_slot, aplus_raw,
        )
        bminus = _conditional_endpoint_slot!(
            ws, plan, terminal.b, plan.terminal_b_slot, bminus_raw,
        )
        bplus = _conditional_endpoint_slot!(
            ws, plan, terminal.b, plan.terminal_b_slot, bplus_raw,
        )
        H = _pair_rectprob_counted!(
            ws, plan.terminal_copula, aminus, aplus, bminus, bplus,
        )
        G = weight * H
    else
        first = ranks[terminal.a] < ranks[terminal.b] ? terminal.a : terminal.b
        second = first == terminal.a ? terminal.b : terminal.a
        first_slot = first == terminal.a ? plan.terminal_a_slot : plan.terminal_b_slot
        raw_lower, raw_upper = _ordinal_gap(first, ranks, conditioning_values, sampled)
        alpha = _conditional_endpoint_slot!(ws, plan, first, first_slot, raw_lower)
        beta = _conditional_endpoint_slot!(ws, plan, first, first_slot, raw_upper)
        q_terminal = clamp(beta - alpha, 0.0, 1.0)
        iszero(q_terminal) && return 0.0, q1
        r = clamp(
            alpha + clamp(Float64(z[end]), 0.0, 1.0) * q_terminal,
            alpha, beta,
        )
        raw_first = _invert_reduced_state!(ws, plan, first_slot, r)
        conditioning_values[first] = _set_raw_value!(ws, plan, first, raw_first)
        second_given_first = plan.states[
            _reduced_state_key(second, vcat(terminal.conditioning, first))
        ]
        upper_conditional = _conditional_endpoint_slot!(
            ws, plan, second, second_given_first, raw_upper,
        )
        lower_conditional = _conditional_endpoint_slot!(
            ws, plan, second, second_given_first, raw_first,
        )
        conditional_probability = clamp(
            upper_conditional - lower_conditional, 0.0, 1.0,
        )
        G = weight * q_terminal * conditional_probability
    end

    tolerance = 4096eps(Float64) * max(1.0, q1)
    (-tolerance <= G <= q1 + tolerance) || throw(DomainError(
        G, "ordinal reduced integrand violates 0 <= G <= q1",
    ))
    return clamp(G, 0.0, q1), q1
end

function _ordinal_probability_rqmc(
    vc::AbstractVineCopula{p},
    pattern::AbstractVector{<:Integer};
    N::Integer=1 << 12,
    R::Integer=8,
    seed::Integer=0x4f52444e,
    sampling_order=nothing,
    randomized::Bool=true,
) where {p}
    Ni, Ri = Int(N), Int(R)
    Ri >= 1 || throw(ArgumentError("R must be positive"))
    ranks = _ordinal_ranks(pattern, p)
    plan = _compile_reduced_probability(vc)
    order = sampling_order === nothing ? copy(plan.default_order) :
        _validate_reduced_order(plan, sampling_order)
    slots = _reduced_sampling_slots(plan, order)
    adjacent = abs(ranks[plan.terminal.a] - ranks[plan.terminal.b]) == 1
    dimension = p - 2 + Int(adjacent)
    replicates = zeros(Float64, Ri)
    ws = _reduced_workspace(plan)
    values = fill(NaN, p)
    sampled = falses(p)
    for rep in 1:Ri
        Z = _reduced_owen_points(
            dimension, Ni, UInt64(seed) + UInt64(rep - 1);
            randomized=randomized,
        )
        total = 0.0
        @inbounds for col in axes(Z, 2)
            G, _ = _ordinal_probability_integrand!(
                ws, plan, pattern, ranks, view(Z, :, col), order, slots,
                values, sampled,
            )
            total += G
        end
        replicates[rep] = total / Ni
    end
    estimate = sum(replicates) / Ri
    stderr = Ri > 1 ? sqrt(sum(abs2, replicates .- estimate) / (Ri - 1)) /
        sqrt(Ri) : NaN
    return _ReducedRQMCResult(
        estimate, stderr, replicates, Ni, Ri, dimension, order,
        ws.counts.hfunc, ws.counts.hinv, ws.counts.pair_cdf,
    )
end

function _rectprob_reduced_rqmc_with_order(
    vc::AbstractVineCopula{p},
    lower::AbstractVector{<:Real},
    upper::AbstractVector{<:Real},
    sampling_order::AbstractVector{<:Integer};
    N::Integer=1 << 12,
    R::Integer=8,
    seed::Integer=0x5245514d,
    randomized::Bool=true,
) where {p}
    Ni, Ri = Int(N), Int(R)
    lo, hi = _check_reduced_bounds(p, lower, upper)
    plan = _compile_reduced_probability(vc)
    order = _validate_reduced_order(plan, sampling_order)
    slots = _reduced_sampling_slots(plan, order)
    replicates = zeros(Float64, Ri)
    ws = _reduced_workspace(plan)
    values = fill(NaN, p)
    for rep in 1:Ri
        Z = _reduced_owen_points(
            p - 2, Ni, UInt64(seed) + UInt64(rep - 1); randomized=randomized,
        )
        total = 0.0
        @inbounds for col in axes(Z, 2)
            G, _ = _reduced_probability_integrand!(
                ws, plan, lo, hi, view(Z, :, col), order, slots, values,
            )
            total += G
        end
        replicates[rep] = total / Ni
    end
    estimate = sum(replicates) / Ri
    stderr = Ri > 1 ?
        sqrt(sum(abs2, replicates .- estimate) / (Ri - 1)) / sqrt(Ri) : NaN
    return _ReducedRQMCResult(
        estimate, stderr, replicates, Ni, Ri, p - 2, order,
        ws.counts.hfunc, ws.counts.hinv, ws.counts.pair_cdf,
    )
end

function _select_reduced_order_pilot(
    vc::AbstractVineCopula{p},
    lower::AbstractVector{<:Real},
    upper::AbstractVector{<:Real};
    M::Integer=1 << 12,
    seed::Integer=0x50494c54,
    candidate_orders=nothing,
) where {p}
    lo, hi = _check_reduced_bounds(p, lower, upper)
    plan = _compile_reduced_probability(vc)
    orders = candidate_orders === nothing ? _reduced_admissible_orders(plan) :
        [_validate_reduced_order(plan, order) for order in candidate_orders]
    isempty(orders) && throw(ArgumentError("no admissible candidate orders"))
    rng = Random.MersenneTwister(seed)
    Z = rand(rng, p - 2, Int(M))
    variances = Vector{Float64}(undef, length(orders))
    means = similar(variances)
    ws = _reduced_workspace(plan)
    values = fill(NaN, p)
    for (index, order) in enumerate(orders)
        slots = _reduced_sampling_slots(plan, order)
        samples = Vector{Float64}(undef, Int(M))
        @inbounds for col in axes(Z, 2)
            samples[col], _ = _reduced_probability_integrand!(
                ws, plan, lo, hi, view(Z, :, col), order, slots, values,
            )
        end
        means[index] = sum(samples) / length(samples)
        variances[index] = length(samples) > 1 ?
            sum(abs2, samples .- means[index]) / (length(samples) - 1) : 0.0
    end
    selected = argmin(variances)
    return (;
        order=orders[selected], selected_index=selected, orders,
        variances, means, M=Int(M), seed=Int(seed),
    )
end

function _rectprob_reduced_pilot_rqmc(
    vc::AbstractVineCopula,
    lower::AbstractVector{<:Real},
    upper::AbstractVector{<:Real};
    M::Integer=1 << 12,
    N::Integer=1 << 12,
    R::Integer=8,
    pilot_seed::Integer=0x50494c54,
    final_seed::Integer=0x46494e4c,
    candidate_orders=nothing,
)
    pilot_seed == final_seed && throw(ArgumentError(
        "pilot and final seeds must differ",
    ))
    pilot = _select_reduced_order_pilot(
        vc, lower, upper; M, seed=pilot_seed, candidate_orders,
    )
    result = _rectprob_reduced_rqmc_with_order(
        vc, lower, upper, pilot.order; N, R, seed=final_seed,
    )
    return (; result, pilot)
end
