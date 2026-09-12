using AdditionalDistributions
using Copulas
using Distributions
using LinearAlgebra
using Printf
using Random
using Statistics
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
const N = 4096
const REPS = 30
const BASE_SEED = 271_828
const AD_BASE_SEED = 729_828
const AD_SHIFTS = 10
const AD_ABSEPS = 1e-6
const AD_RELEPS = 1e-6

"""Recover the Gaussian correlation matrix from D-vine partial correlations."""
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

function load_references()
    rows = read_table(joinpath(FROZEN_DIR, "references.csv"))
    return Dict(row["case_id"] => row for row in rows)
end

function load_frozen_raw()
    rows = read_table(joinpath(FROZEN_DIR, "raw_results.csv"))
    lookup = Dict{Tuple{String,Int,String},Float64}()
    for row in rows
        _parse_int(row, "N") == N || continue
        row["case_id"] in SELECTED_CASES || continue
        method = row["method"]
        method in ("current_indicator", "reduced_conditional") || continue
        lookup[(row["case_id"], _parse_int(row, "rep"), method)] = _parse_float(row, "estimate")
    end
    return lookup
end

function selected_cases()
    cfg = (dims=[4, 5], case_set="core")
    cases = filter(case -> case.id in SELECTED_CASES, publication_cases(cfg))
    sort!(cases; by=case -> case.id)
    length(cases) == length(SELECTED_CASES) || error("not all selected cases were found")
    return cases
end

function method_label(method)
    method == "current_indicator" && return "Indicator"
    method == "reduced_conditional" && return "Reduced"
    method == "additionaldistributions" && return "Gaussian specialized"
    error("unknown method $method")
end

function event_label(case)
    case.point == "rectangle" && return "interior"
    return "CDF"
end

function dependence_label(case)
    occursin("strong", case.id) && return "strong"
    return "moderate"
end

function run_comparison()
    refs = load_references()
    frozen = load_frozen_raw()
    cases = selected_cases()
    raw = NamedTuple[]
    max_parity_error = 0.0
    max_frozen_difference = 0.0

    for case in cases
        reference = _parse_float(refs[case.id], "estimate")
        R = gaussian_dvine_correlation(case.vine)
        gaussian_copula = GaussianCopula(R)
        probe = collect(range(0.19, 0.81; length=length(case.vine)))
        parity_error = abs(logpdf(case.vine, probe) - logpdf(gaussian_copula, probe))
        max_parity_error = max(max_parity_error, parity_error)
        parity_error <= 2e-11 || error("Gaussian density parity failed for $(case.id): $parity_error")

        lower_z = normal_limit.(case.lower)
        upper_z = normal_limit.(case.upper)
        mvn = MvNormal(zeros(length(case.vine)), Symmetric(R))
        runner = prepare_reduced_runner(case)

        # Exclude JIT, Reduced-plan construction, and first CBC-lattice construction.
        current_estimate(case.vine, case.lower, case.upper, 64, BASE_SEED)
        reduced_compiled_estimate!(runner, 64, BASE_SEED)
        AdditionalDistributions.cdf_result(
            mvn, lower_z, upper_z;
            m=N, abseps=AD_ABSEPS, releps=AD_RELEPS,
            nshifts=AD_SHIFTS, rng=MersenneTwister(AD_BASE_SEED),
        )

        for rep in 1:REPS
            order = mod1(rep, 3)
            methods = order == 1 ? ("current_indicator", "reduced_conditional", "additionaldistributions") :
                      order == 2 ? ("reduced_conditional", "additionaldistributions", "current_indicator") :
                                   ("additionaldistributions", "current_indicator", "reduced_conditional")

            for method in methods
                internal_error = NaN
                inform = -1
                algorithm = ""
                seed = method == "additionaldistributions" ? AD_BASE_SEED + rep - 1 : BASE_SEED + rep - 1
                if method == "current_indicator"
                    timed = @timed current_estimate(case.vine, case.lower, case.upper, N, seed)
                    estimate = timed.value
                elseif method == "reduced_conditional"
                    timed = @timed reduced_compiled_estimate!(runner, N, seed)
                    estimate = timed.value.estimate
                else
                    timed = @timed AdditionalDistributions.cdf_result(
                        mvn, lower_z, upper_z;
                        m=N, abseps=AD_ABSEPS, releps=AD_RELEPS,
                        nshifts=AD_SHIFTS, rng=MersenneTwister(seed),
                    )
                    result = timed.value
                    estimate = result.value
                    internal_error = result.error
                    inform = result.inform
                    algorithm = string(result.algorithm)
                end

                if method != "additionaldistributions"
                    key = (case.id, rep, method)
                    haskey(frozen, key) || error("missing frozen result for $key")
                    max_frozen_difference = max(max_frozen_difference, abs(estimate - frozen[key]))
                end

                push!(raw, (;
                    case_id=case.id,
                    d=length(case.vine),
                    event=event_label(case),
                    dependence=dependence_label(case),
                    method,
                    N_or_maxpts=N,
                    replicate=rep,
                    seed,
                    estimate,
                    reference,
                    abs_error=abs(estimate - reference),
                    runtime_s=timed.time,
                    internal_error,
                    inform,
                    algorithm,
                    nshifts=method == "additionaldistributions" ? AD_SHIFTS : 1,
                ))
            end
        end
    end

    max_frozen_difference == 0.0 || error("frozen Reduced/Indicator estimates changed: $max_frozen_difference")
    write_rows(joinpath(OUT_DIR, "raw_results.csv"), raw)
    return raw, max_parity_error, max_frozen_difference
end

function summarize(raw)
    groups = Dict{Tuple{String,String},Vector{typeof(first(raw))}}()
    for row in raw
        push!(get!(groups, (row.case_id, row.method), typeof(first(raw))[]), row)
    end
    summary = NamedTuple[]
    method_rank = Dict(
        "current_indicator" => 1,
        "reduced_conditional" => 2,
        "additionaldistributions" => 3,
    )
    ordered_keys = sort(collect(keys(groups)); by = key -> (
        first(groups[key]).d,
        key[1],
        method_rank[key[2]],
    ))
    for key in ordered_keys
        rows = groups[key]
        firstrow = first(rows)
        errors = getproperty.(rows, :abs_error)
        times = getproperty.(rows, :runtime_s)
        internal = filter(isfinite, getproperty.(rows, :internal_error))
        push!(summary, (;
            case_id=firstrow.case_id,
            d=firstrow.d,
            event=firstrow.event,
            dependence=firstrow.dependence,
            method=firstrow.method,
            N_or_maxpts=firstrow.N_or_maxpts,
            replicates=length(rows),
            reference=firstrow.reference,
            rmse=sqrt(mean(abs2, getproperty.(rows, :estimate) .- firstrow.reference)),
            median_abs_error=median(errors),
            median_runtime_s=median(times),
            median_internal_error=isempty(internal) ? NaN : median(internal),
            tolerance=firstrow.method == "additionaldistributions" ? AD_ABSEPS : NaN,
            successful_statuses=count(==(0), getproperty.(rows, :inform)),
        ))
    end
    write_rows(joinpath(OUT_DIR, "summary.csv"), summary)
    return summary
end

raw, parity, frozen_difference = run_comparison()
summary = summarize(raw)

println("Gaussian specialized comparison complete")
println("  cases                         = ", length(SELECTED_CASES))
println("  rows                          = ", length(raw))
println("  maximum density parity error = ", parity)
println("  maximum frozen estimate diff = ", frozen_difference)
println("  AdditionalDistributions      = ", pkgversion(AdditionalDistributions))
println("  raw                           = ", joinpath(OUT_DIR, "raw_results.csv"))
println("  summary                       = ", joinpath(OUT_DIR, "summary.csv"))
