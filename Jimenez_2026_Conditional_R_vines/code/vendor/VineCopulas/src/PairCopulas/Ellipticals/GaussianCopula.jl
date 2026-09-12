# ---------------------------------------------------------------------
# Gaussian pair-copula fast path
# ---------------------------------------------------------------------

@inline function _gaussian_pair_inputs(C::Copulas.GaussianCopula{2}, u::Real, v::Real)
    ρ = C.Σ[1, 2]
    z1 = Distributions.quantile(_STD_NORMAL, _clp(u))
    z2 = Distributions.quantile(_STD_NORMAL, _clp(v))
    ρ2 = ρ * ρ
    den = one(ρ2) - ρ2
    return ρ, z1, z2, den
end

@inline function _gaussian_logpdf_from_z(ρ::Real, z1::Real, z2::Real, den::Real)
    ρ2 = ρ * ρ
    return -0.5 * log(den) +
           (2 * ρ * z1 * z2 - ρ2 * (z1 * z1 + z2 * z2)) / (2 * den)
end

@inline _gaussian_h_from_z(ρ::Real, target_z::Real, base_z::Real, den::Real) =
    Distributions.cdf(_STD_NORMAL, (target_z - ρ * base_z) / sqrt(den))

@inline function _pair_logpdf(
    C::Copulas.GaussianCopula{2},
    u::Real,
    v::Real,
    buf::Vector{Float64},
)
    ρ, z1, z2, den = _gaussian_pair_inputs(C, u, v)
    return _gaussian_logpdf_from_z(ρ, z1, z2, den)
end

# Copulas.jl currently routes even the bivariate Gaussian CDF through a
# randomized multivariate-normal QMC routine. The exact deterministic
# bivariate routine already provided by StatsFuns is essential when this CDF
# is evaluated inside an outer QMC integrand.
@inline function _pair_cdf(
    C::Copulas.GaussianCopula{2},
    u::Real,
    v::Real,
    buf::Vector{Float64},
)
    uu, vv = clamp(Float64(u), 0.0, 1.0), clamp(Float64(v), 0.0, 1.0)
    (iszero(uu) || iszero(vv)) && return 0.0
    isone(uu) && return vv
    isone(vv) && return uu
    z1 = Distributions.quantile(_STD_NORMAL, uu)
    z2 = Distributions.quantile(_STD_NORMAL, vv)
    return clamp(StatsFuns.bvncdf(z1, z2, Float64(C.Σ[1, 2])), 0.0, 1.0)
end

@inline _pair_cdf_safety(::Copulas.GaussianCopula{2}) = :deterministic_analytic

@inline function _pair_hfuncs(C::Copulas.GaussianCopula{2}, u::Real, v::Real)
    ρ, z1, z2, den = _gaussian_pair_inputs(C, u, v)
    h1 = _clp(_gaussian_h_from_z(ρ, z1, z2, den))
    h2 = _clp(_gaussian_h_from_z(ρ, z2, z1, den))
    return h1, h2
end

@inline function _pair_step(
    C::Copulas.GaussianCopula{2},
    u::Real,
    v::Real,
    buf::Vector{Float64},
)
    ρ, z1, z2, den = _gaussian_pair_inputs(C, u, v)
    logc = _gaussian_logpdf_from_z(ρ, z1, z2, den)
    h1 = _clp(_gaussian_h_from_z(ρ, z1, z2, den))
    h2 = _clp(_gaussian_h_from_z(ρ, z2, z1, den))
    return logc, h1, h2
end

@inline function _pair_logpdf_h1(
    C::Copulas.GaussianCopula{2},
    u::Real,
    v::Real,
    buf::Vector{Float64},
)
    ρ, z1, z2, den = _gaussian_pair_inputs(C, u, v)
    logc = _gaussian_logpdf_from_z(ρ, z1, z2, den)
    h1 = _clp(_gaussian_h_from_z(ρ, z1, z2, den))
    return logc, h1
end

@inline function _pair_logpdf_h2(
    C::Copulas.GaussianCopula{2},
    u::Real,
    v::Real,
    buf::Vector{Float64},
)
    ρ, z1, z2, den = _gaussian_pair_inputs(C, u, v)
    logc = _gaussian_logpdf_from_z(ρ, z1, z2, den)
    h2 = _clp(_gaussian_h_from_z(ρ, z2, z1, den))
    return logc, h2
end

@inline function _gaussian_hfunc(C::Copulas.GaussianCopula{2}, target::Real, base::Real)
    ρ = C.Σ[1, 2]
    zt = Distributions.quantile(_STD_NORMAL, _clp(target))
    zb = Distributions.quantile(_STD_NORMAL, _clp(base))
    den = one(ρ) - ρ * ρ
    return _gaussian_h_from_z(ρ, zt, zb, den)
end

@inline function _gaussian_hinv(C::Copulas.GaussianCopula{2}, q::Real, base::Real)
    ρ = C.Σ[1, 2]
    zb = Distributions.quantile(_STD_NORMAL, _clp(base))
    zq = Distributions.quantile(_STD_NORMAL, _clp(q))
    return Distributions.cdf(_STD_NORMAL, ρ * zb + sqrt(1 - ρ^2) * zq)
end

@inline hfunc1(C::Copulas.GaussianCopula{2}, u::Real, v::Real) = _clp(_gaussian_hfunc(C, u, v))
@inline hfunc2(C::Copulas.GaussianCopula{2}, u::Real, v::Real) = _clp(_gaussian_hfunc(C, v, u))

@inline hfunc1(C::Copulas.GaussianCopula{2}, uv::Tuple{<:Real,<:Real}) = hfunc1(C, uv[1], uv[2])

@inline hfunc2(C::Copulas.GaussianCopula{2}, uv::Tuple{<:Real,<:Real}) = hfunc2(C, uv[1], uv[2])

@inline hinv1(C::Copulas.GaussianCopula{2}, q::Real, v::Real) = _clp(_gaussian_hinv(C, q, v))
@inline hinv2(C::Copulas.GaussianCopula{2}, q::Real, u::Real) = _clp(_gaussian_hinv(C, q, u))
