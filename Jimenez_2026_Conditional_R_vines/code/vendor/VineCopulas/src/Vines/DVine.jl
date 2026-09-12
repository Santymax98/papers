# D-vine pair-copula constructions.
# Convention for edges[k][i]: C_{left, right | inner chain}, with
# left = order[i], right = order[i+k] and conditioning set order[i+1:i+k-1].
# Copula coordinates are (left, right).

"""Structure-only description of a D-vine path order and active truncation level."""
struct DVineStructure{p,q} <: AbstractVineStructure{p}
    order::NTuple{p,Int}

    function DVineStructure{p,q}(order::NTuple{p,Int}) where {p,q}
        1 <= q <= p - 1 || throw(ArgumentError("trunc must be in 1:$(p-1)"))
        sort(collect(order)) == collect(1:p) ||
            throw(ArgumentError("order must be a permutation of 1:$p"))
        return new{p,q}(order)
    end
end

function DVineStructure(order::_OrderInput; trunc::Int=length(order)-1)
    p = _check_order(order)
    1 <= trunc <= p-1 || throw(ArgumentError("trunc must be in 1:$(p-1)"))
    return DVineStructure{p,trunc}(Tuple(Int.(order)))
end

"""
    DVineCopula(order, edges; trunc=length(order)-1)
    DVineCopula(structure::DVineStructure, edges)

Construct a drawable/path vine copula from a variable `order` and a triangular
collection of bivariate pair-copulas. The entry `edges[k][i]` represents the
pair-copula between `order[i]` and `order[i+k]`, conditional on the variables
between them in the D-vine path.

# Example

```julia
C12 = GaussianCopula([1.0 0.5; 0.5 1.0])
C23 = ClaytonCopula(2, 2.0)
C13_2 = FrankCopula(2, 3.0)
dv = DVineCopula([1, 2, 3], [[C12, C23], [C13_2]])
```
"""
struct DVineCopula{p,q,E} <: AbstractVineCopula{p}
    order::NTuple{p,Int}
    edges::E
    trunc::Int
end

function DVineCopula(; order, paircopulas, trunc = length(order) - 1)
    ord = collect(Int, order)
    pcs = [collect(level) for level in paircopulas]
    return DVineCopula(ord, pcs; trunc = trunc)
end

function DVineCopula(order::_OrderInput, edges; trunc::Int=length(order)-1)
    p = _check_order(order)
    1 <= trunc <= p-1 || throw(ArgumentError("trunc must be in 1:$(p-1)"))
    E = _normalize_edges(edges, p, trunc)
    return DVineCopula{p,trunc,typeof(E)}(Tuple(Int.(order)), E, trunc)
end

function DVineCopula(structure::DVineStructure{p,q}, edges) where {p,q}
    E = _normalize_edges(edges, p, q)
    return DVineCopula{p,q,typeof(E)}(structure.order, E, q)
end

function DVineCopula(edges; order=nothing, trunc::Int=length(edges))
    p = length(edges) + 1
    order === nothing && (order = collect(1:p))
    return DVineCopula(order, edges; trunc=trunc)
end

"""
    order(vine)

Return the variable order used by a vine copula.
"""
order(vc::DVineCopula) = vc.order
order(st::DVineStructure) = st.order

"""Return a `DVineStructure` describing the D-vine path order and truncation."""
structure(vc::DVineCopula{p,q}) where {p,q} = DVineStructure{p,q}(vc.order)

"""
    edges(vine)

Return the triangular array of bivariate pair-copulas used by a vine copula.
Tree `k` is stored in `edges(vine)[k]`.
"""
edges(vc::DVineCopula) = vc.edges

"""
    truncation(vine)

Return the number of active trees in the vine. A full `p`-dimensional vine has
truncation level `p - 1`.
"""
truncation(vc::DVineCopula) = vc.trunc
truncation(::DVineStructure{p,q}) where {p,q} = q

function truncate(st::DVineStructure{p}, level::Integer) where {p}
    q = _check_truncate_level(level, p, truncation(st))
    return DVineStructure{p,q}(st.order)
end

function truncate(vc::DVineCopula{p}, level::Integer) where {p}
    st = truncate(structure(vc), level)
    return DVineCopula(st, vc.edges[1:truncation(st)])
end

Base.show(io::IO, vc::DVineCopula{p}) where {p} = print(io, "DVineCopula(p=$p, trunc=$(vc.trunc))")
Base.show(io::IO, st::DVineStructure{p}) where {p} = print(io, "DVineStructure(p=$p, trunc=$(truncation(st)))")

function _logpdf_internal(vc::DVineCopula{p}, u::AbstractVector{<:Real}) where {p}
    _check_vector_dim(p, u)
    return _logpdf_internal(vc, reshape(u, p, 1))[1]
end

function _logpdf_internal(vc::DVineCopula{p}, U::AbstractMatrix{<:Real}) where {p}
    X = _as_pxn(p, U)
    n = size(X,2)
    L = Matrix{Float64}(undef, p, n)
    R = Matrix{Float64}(undef, p, n)
    @inbounds for j in 1:p
        @views L[j,:] .= _clp.(X[vc.order[j],:])
        @views R[j,:] .= L[j,:]
    end
    ll = zeros(Float64, n)
    buf = Vector{Float64}(undef, 2)
    @inbounds for k in 1:vc.trunc
        # Pair (i, i+k | i+1:i+k-1) uses L[i] and R[i+k].  Within one
        # tree every written L[i] / R[i+k] slot is unique, so conditionals can
        # replace the consumed states in place without copying full matrices.
        propagate = k < vc.trunc
        for i in 1:(p-k)
            C = vc.edges[k][i]
            left = @view L[i,:]
            right = @view R[i+k,:]
            if propagate
                # Dispatch on a heterogeneous pair container happens once at
                # this function barrier. Both input rows may alias the outputs.
                _pair_step_add!(ll, left, right, C, left, right, buf)
            else
                _pair_logpdf_add!(ll, C, left, right, buf)
            end
        end
    end
    return ll
end

function _dvine_left_conditionals!(
    Lwork::Matrix{Float64},
    Rwork::Vector{Float64},
    vc::DVineCopula{p},
    X::AbstractMatrix{<:Real},
    i::Int,
) where {p}
    # Lwork[m,:] = u_{m | m+1:(i-1)} is needed only in the active truncation
    # window. A single Rwork vector carries the right conditional down each
    # recurrence, replacing the previous vector-of-vectors construction.
    n = size(X, 2)
    size(Lwork) == (p, n) || throw(DimensionMismatch("Lwork has incompatible size"))
    length(Rwork) == n || throw(DimensionMismatch("Rwork has incompatible length"))
    first = max(1, i - vc.trunc)

    @inbounds for m in first:(i-1), col in 1:n
        Lwork[m,col] = X[m,col]
    end

    @inbounds for t in (first+1):(i-1)
        for col in 1:n
            Rwork[col] = X[t,col]
        end
        for m in (t-1):-1:first
            C = vc.edges[t-m][m]
            for col in 1:n
                uL = _clp(Lwork[m,col])
                uR = _clp(Rwork[col])
                h1, h2 = _pair_hfuncs(C, uL, uR)
                Lwork[m,col] = h1
                Rwork[col] = h2
            end
        end
    end
    return first
end

function _rosenblatt_internal!(out::AbstractMatrix{<:Real}, vc::DVineCopula{p}, U::AbstractMatrix{<:Real}) where {p}
    Ux = _as_pxn(p, U)
    n = size(Ux,2)
    X = Matrix{Float64}(undef, p, n)
    @inbounds for j in 1:p
        @views X[j,:] .= _clp.(Ux[vc.order[j],:])
    end
    Z = Matrix{Float64}(undef, p, n)
    Lwork = Matrix{Float64}(undef, p, n)
    Rwork = Vector{Float64}(undef, n)
    @inbounds Z[1,:] .= X[1,:]
    @inbounds for i in 2:p
        first = _dvine_left_conditionals!(Lwork, Rwork, vc, X, i)
        @views Z[i,:] .= X[i,:]
        # z_i = F_{i | 1:(i-1)}. Apply hfunc2 from nearest to farthest left.
        for m in (i-1):-1:first
            C = vc.edges[i-m][m]
            for col in 1:n
                Z[i,col] = hfunc2(C, Lwork[m,col], Z[i,col])
            end
        end
    end
    invord = _invperm_tuple(vc.order)
    @inbounds for label in 1:p
        @views out[label,:] .= Z[invord[label],:]
    end
    return out
end

function _inverse_rosenblatt_internal!(out::AbstractMatrix{<:Real}, vc::DVineCopula{p}, Z::AbstractMatrix{<:Real}) where {p}
    Zx = _as_pxn(p, Z)
    n = size(Zx,2)
    W = Matrix{Float64}(undef, p, n)
    @inbounds for j in 1:p
        @views W[j,:] .= _clp.(Zx[vc.order[j],:])
    end
    X = Matrix{Float64}(undef, p, n)
    Lwork = Matrix{Float64}(undef, p, n)
    Rwork = Vector{Float64}(undef, n)
    @inbounds X[1,:] .= W[1,:]
    @inbounds for i in 2:p
        first = _dvine_left_conditionals!(Lwork, Rwork, vc, X, i)
        @views X[i,:] .= W[i,:]
        # Invert from farthest conditioned pair to nearest.
        for m in first:(i-1)
            C = vc.edges[i-m][m]
            for col in 1:n
                X[i,col] = hinv2(C, X[i,col], Lwork[m,col])
            end
        end
    end
    invord = _invperm_tuple(vc.order)
    @inbounds for label in 1:p
        @views out[label,:] .= X[invord[label],:]
    end
    return out
end

# ---------------------------------------------------------------------
# Experimental reduced CDF integrals for the full four-dimensional
# D-vine with path order 1-2-3-4. These functions are deliberately private:
# they test a numerical representation without changing the public CDF API.
# ---------------------------------------------------------------------

@inline function _check_dvine4_reduced(vc::DVineCopula{4}, u::AbstractVector{<:Real})
    vc.order == (1, 2, 3, 4) ||
        throw(ArgumentError("the experimental reduced CDF requires D-vine order (1, 2, 3, 4)"))
    vc.trunc == 3 ||
        throw(ArgumentError("the experimental reduced CDF requires a full four-dimensional D-vine"))
    _check_vector_dim(4, u)
    any(isnan, u) && throw(DomainError(u, "CDF coordinates cannot be NaN"))
    return ntuple(i -> clamp(Float64(u[i]), 0.0, 1.0), 4)
end

@inline function _dvine4_reduced_edges(vc::DVineCopula{4})
    C12, C23, C34 = vc.edges[1]
    C13_2, C24_3 = vc.edges[2]
    C14_23 = vc.edges[3][1]
    return C12, C23, C34, C13_2, C24_3, C14_23
end

@inline function _dvine4_reduced_affine_value(
    C12::C12T,
    C23::C23T,
    C34::C34T,
    C13_2::C13T,
    C24_3::C24T,
    C14_23::C14T,
    u::NTuple{4,Float64},
    z1::Real,
    z2::Real,
    pdfbuf::Vector{Float64},
    cdfbuf::Vector{Float64},
) where {C12T<:PairCopula,C23T<:PairCopula,C34T<:PairCopula,
         C13T<:PairCopula,C24T<:PairCopula,C14T<:PairCopula}
    u1, u2, u3, u4 = u
    v2, v3 = u2 * z1, u3 * z2

    r1_2 = hfunc1(C12, u1, v2)
    r2_3, r3_2 = _pair_hfuncs(C23, v2, v3)
    A = hfunc1(C13_2, r1_2, r3_2)
    r4_3 = hfunc2(C34, v3, u4)
    B = hfunc2(C24_3, r2_3, r4_3)

    c23 = exp(_pair_logpdf(C23, v2, v3, pdfbuf))
    return u2 * u3 * c23 * _pair_cdf(C14_23, A, B, cdfbuf)
end

@inline function _dvine4_reduced_conditional_value(
    C12::C12T,
    C23::C23T,
    C34::C34T,
    C13_2::C13T,
    C24_3::C24T,
    C14_23::C14T,
    u::NTuple{4,Float64},
    z1::Real,
    z2::Real,
    cdfbuf::Vector{Float64},
) where {C12T<:PairCopula,C23T<:PairCopula,C34T<:PairCopula,
         C13T<:PairCopula,C24T<:PairCopula,C14T<:PairCopula}
    u1, u2, u3, u4 = u
    v2 = u2 * z1
    q3 = hfunc2(C23, v2, u3)
    r3_2 = z2 * q3
    v3 = hinv2(C23, r3_2, v2)

    r1_2 = hfunc1(C12, u1, v2)
    A = hfunc1(C13_2, r1_2, r3_2)
    r2_3 = hfunc1(C23, v2, v3)
    r4_3 = hfunc2(C34, v3, u4)
    B = hfunc2(C24_3, r2_3, r4_3)

    return u2 * q3 * _pair_cdf(C14_23, A, B, cdfbuf)
end

"""Private scalar integrand used to verify the experimental affine representation."""
function _dvine4_reduced_affine_integrand(
    vc::DVineCopula{4},
    u::AbstractVector{<:Real},
    z::AbstractVector{<:Real},
)
    length(z) == 2 || throw(ArgumentError("z must have length 2"))
    uc = _check_dvine4_reduced(vc, u)
    E = _dvine4_reduced_edges(vc)
    return _dvine4_reduced_affine_value(E..., uc, z[1], z[2], zeros(2), zeros(2))
end

"""Private scalar integrand used to verify the experimental conditional representation."""
function _dvine4_reduced_conditional_integrand(
    vc::DVineCopula{4},
    u::AbstractVector{<:Real},
    z::AbstractVector{<:Real},
)
    length(z) == 2 || throw(ArgumentError("z must have length 2"))
    uc = _check_dvine4_reduced(vc, u)
    E = _dvine4_reduced_edges(vc)
    return _dvine4_reduced_conditional_value(E..., uc, z[1], z[2], zeros(2))
end

@inline function _dvine4_reduced_boundary(vc::DVineCopula{4}, u::NTuple{4,Float64})
    any(iszero, u) && return 0.0
    all(isone, u) && return 1.0
    if isone(u[1]) && isone(u[4])
        return Float64(_pair_cdf(vc.edges[1][2], u[2], u[3], zeros(2)))
    end
    return nothing
end


"""
    _cdf_dvine4_reduced_affine(vc, u; N=_CDF_NSAMPLES[], randomized=true)

Experimental QMC approximation of a four-dimensional D-vine CDF using the
two-dimensional affine integral over `(v2, v3)`. This is intentionally private.
"""
function _cdf_dvine4_reduced_affine(
    vc::DVineCopula{4},
    u::AbstractVector{<:Real};
    N::Integer=_CDF_NSAMPLES[],
    randomized::Bool=true,
)
    N >= 1 || throw(ArgumentError("N must be positive"))
    uc = _check_dvine4_reduced(vc, u)
    boundary = _dvine4_reduced_boundary(vc, uc)
    boundary === nothing || return boundary
    Z = _qmc_points(2, Int(N); randomized=randomized)
    E = _dvine4_reduced_edges(vc)
    return _cdf_dvine4_reduced_affine_qmc(E..., uc, Z)
end

function _cdf_dvine4_reduced_affine_qmc(
    C12::C12T,
    C23::C23T,
    C34::C34T,
    C13_2::C13T,
    C24_3::C24T,
    C14_23::C14T,
    u::NTuple{4,Float64},
    Z::AbstractMatrix{<:Real},
) where {C12T<:PairCopula,C23T<:PairCopula,C34T<:PairCopula,
         C13T<:PairCopula,C24T<:PairCopula,C14T<:PairCopula}
    pdfbuf, cdfbuf = zeros(2), zeros(2)
    total = 0.0
    @inbounds for j in axes(Z, 2)
        total += _dvine4_reduced_affine_value(
            C12, C23, C34, C13_2, C24_3, C14_23,
            u, Z[1, j], Z[2, j], pdfbuf, cdfbuf,
        )
    end
    return total / size(Z, 2)
end

"""
    _cdf_dvine4_reduced_conditional(vc, u; N=_CDF_NSAMPLES[], randomized=true)

Experimental QMC approximation of a four-dimensional D-vine CDF using the
two-dimensional `2 -> 3` conditional transform. This is intentionally private.
"""
function _cdf_dvine4_reduced_conditional(
    vc::DVineCopula{4},
    u::AbstractVector{<:Real};
    N::Integer=_CDF_NSAMPLES[],
    randomized::Bool=true,
)
    N >= 1 || throw(ArgumentError("N must be positive"))
    uc = _check_dvine4_reduced(vc, u)
    boundary = _dvine4_reduced_boundary(vc, uc)
    boundary === nothing || return boundary
    Z = _qmc_points(2, Int(N); randomized=randomized)
    E = _dvine4_reduced_edges(vc)
    return _cdf_dvine4_reduced_conditional_qmc(E..., uc, Z)
end

function _cdf_dvine4_reduced_conditional_qmc(
    C12::C12T,
    C23::C23T,
    C34::C34T,
    C13_2::C13T,
    C24_3::C24T,
    C14_23::C14T,
    u::NTuple{4,Float64},
    Z::AbstractMatrix{<:Real},
) where {C12T<:PairCopula,C23T<:PairCopula,C34T<:PairCopula,
         C13T<:PairCopula,C24T<:PairCopula,C14T<:PairCopula}
    cdfbuf = zeros(2)
    total = 0.0
    @inbounds for j in axes(Z, 2)
        total += _dvine4_reduced_conditional_value(
            C12, C23, C34, C13_2, C24_3, C14_23,
            u, Z[1, j], Z[2, j], cdfbuf,
        )
    end
    return total / size(Z, 2)
end

function _dvine_edge_description(vc::DVineCopula{p}, k::Int, i::Int) where {p}
    left = vc.order[i]
    right = vc.order[i+k]
    D = Tuple(vc.order[i+1:i+k-1])
    return VineEdge((left, right), D, vc.edges[k][i], k, i)
end

function vine_edges(vc::DVineCopula)
    out = VineEdge[]
    for k in 1:vc.trunc, i in 1:length(vc.edges[k])
        push!(out, _dvine_edge_description(vc, k, i))
    end
    return out
end
