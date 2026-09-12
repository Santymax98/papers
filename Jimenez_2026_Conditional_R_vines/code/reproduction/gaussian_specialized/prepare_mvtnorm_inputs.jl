using Copulas
using Distributions
using LinearAlgebra
using VineCopulas

const HERE = @__DIR__
const PACKAGE_ROOT = normpath(joinpath(HERE, "..", "..", ".."))
const CODE_ROOT = joinpath(PACKAGE_ROOT, "code")
const PUBLICATION_COMMON = joinpath(CODE_ROOT, "benchmarks", "publication", "publication_common.jl")
const FROZEN_DIR = joinpath(CODE_ROOT, "internal_benchmark", "frozen_results")
const OUT_DIR = normpath(get(ENV, "OUT_DIR", joinpath(PACKAGE_ROOT, "output", "gaussian_specialized")))
mkpath(OUT_DIR)

include(PUBLICATION_COMMON)

const SELECTED_CASES = Set([
    "gaussian_moderate_d4",
    "gaussian_strong_d4",
    "gaussian_rectangle_d5",
    "gaussian_strong_d5",
])

function gaussian_dvine_correlation(vine::DVineCopula)
    d = length(vine)
    vine.order == Tuple(1:d) || error("comparison expects D-vine order 1:d")
    vine.trunc == d - 1 || error("comparison expects a full D-vine")
    R = Matrix{Float64}(I, d, d)

    for lag in 1:(d - 1), i in 1:(d - lag)
        j = i + lag
        pair = vine.edges[lag][i]
        pair isa GaussianCopula || error("non-Gaussian pair in selected Gaussian case")
        partial = Float64(pair.Σ[1, 2])
        if lag == 1
            rho = partial
        else
            D = (i + 1):(j - 1)
            RDD = @view R[D, D]
            ri = Vector(@view R[i, D])
            rj = Vector(@view R[j, D])
            xi = RDD \ ri
            xj = RDD \ rj
            projected = dot(ri, xj)
            residual_i = max(0.0, 1.0 - dot(ri, xi))
            residual_j = max(0.0, 1.0 - dot(rj, xj))
            rho = projected + partial * sqrt(residual_i * residual_j)
        end
        R[i, j] = rho
        R[j, i] = rho
    end
    cholesky(Symmetric(R); check=true)
    return R
end

normal_limit(p::Real) = p == 0 ? -Inf : p == 1 ? Inf : quantile(Normal(), p)

references = Dict(row["case_id"] => row for row in read_table(joinpath(FROZEN_DIR, "references.csv")))
cases = filter(case -> case.id in SELECTED_CASES, publication_cases((dims=[4, 5], case_set="core")))
sort!(cases; by=case -> case.id)
length(cases) == length(SELECTED_CASES) || error("not all selected cases were found")

rows = NamedTuple[]
for case in cases
    R = gaussian_dvine_correlation(case.vine)
    probe = collect(range(0.19, 0.81; length=length(case.vine)))
    parity_error = abs(logpdf(case.vine, probe) - logpdf(GaussianCopula(R), probe))
    parity_error <= 2e-11 || error("Gaussian density parity failed for $(case.id): $parity_error")
    push!(rows, (;
        case_id=case.id,
        d=length(case.vine),
        event=case.point == "rectangle" ? "interior" : "CDF",
        dependence=occursin("strong", case.id) ? "strong" : "moderate",
        reference=_parse_float(references[case.id], "estimate"),
        lower_z=_vector_field(normal_limit.(case.lower)),
        upper_z=_vector_field(normal_limit.(case.upper)),
        correlation=_vector_field(vec(R)),
        density_parity_error=parity_error,
    ))
end

output = joinpath(OUT_DIR, "mvtnorm_inputs.csv")
write_rows(output, rows)
println("mvtnorm inputs written to ", output)

