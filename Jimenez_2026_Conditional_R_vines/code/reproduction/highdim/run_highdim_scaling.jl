using Copulas
using Distributions
using LinearAlgebra
using Printf
using QuasiMonteCarlo
using Random
using Statistics
using StatsFuns
using VineCopulas

include(joinpath(@__DIR__, "..", "..", "vendor", "VineCopulas", "benchmarks", "publication", "publication_common.jl"))

const OUTDIR = joinpath(@__DIR__, "rerun")
mkpath(OUTDIR)
const RAW_OUT = joinpath(OUTDIR, "highdim_raw.csv")
const STRESS_OUT = joinpath(OUTDIR, "highdim_stress_raw.csv")
const HALTON_SOBOL_OUT = joinpath(OUTDIR, "highdim_halton_sobol.csv")
const REF_OUT = joinpath(OUTDIR, "highdim_reference_values.csv")
const SUMMARY_OUT = joinpath(OUTDIR, "highdim_summary.csv")
const SKIP_OUT = joinpath(OUTDIR, "highdim_skipped_cells.csv")
const PARITY_OUT = joinpath(OUTDIR, "highdim_density_parity.csv")
const BASE_SEED = parse(Int, get(ENV, "BASE_SEED", "515151"))
const PRIMARY_REPS = parse(Int, get(ENV, "PRIMARY_REPS", "30"))
const STRESS_REPS = parse(Int, get(ENV, "STRESS_REPS", "5"))
const MAX_CELL_SECONDS = parse(Float64, get(ENV, "MAX_CELL_SECONDS", "180"))
const PILOT_ONLY = lowercase(get(ENV, "PILOT_ONLY", "false")) in ("1", "true", "yes")
const FINALIZE_ONLY = lowercase(get(ENV, "FINALIZE_ONLY", "false")) in ("1", "true", "yes")
const SKIP_STRESS_INDICATOR = lowercase(get(ENV, "SKIP_STRESS_INDICATOR", "true")) in ("1", "true", "yes")
const RUN_STRESS = lowercase(get(ENV, "RUN_STRESS", "true")) in ("1", "true", "yes")
const RUN_HALTON_SOBOL = lowercase(get(ENV, "RUN_HALTON_SOBOL", "true")) in ("1", "true", "yes")

csv_escape(x) = begin
    s = string(x)
    (occursin(',', s) || occursin('"', s) || occursin('\n', s)) ?
        "\"" * replace(s, "\"" => "\"\"") * "\"" : s
end

function csv_parse_line(line::AbstractString)
    fields = String[]
    buf = IOBuffer()
    inq = false
    i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if inq
            if c == '"'
                ni = nextind(line, i)
                if ni <= lastindex(line) && line[ni] == '"'
                    write(buf, '"')
                    i = ni
                else
                    inq = false
                end
            else
                write(buf, c)
            end
        elseif c == '"'
            inq = true
        elseif c == ','
            push!(fields, String(take!(buf)))
        else
            write(buf, c)
        end
        i = nextind(line, i)
    end
    push!(fields, String(take!(buf)))
    return fields
end

function read_csv(path)
    isfile(path) || return Dict{String,String}[]
    lines = readlines(path)
    isempty(lines) && return Dict{String,String}[]
    header = csv_parse_line(first(lines))
    rows = Dict{String,String}[]
    for line in Iterators.drop(lines, 1)
        isempty(strip(line)) && continue
        push!(rows, Dict(zip(header, csv_parse_line(line))))
    end
    return rows
end

function append_row(path, row)
    exists = isfile(path)
    open(path, "a") do io
        names = propertynames(row)
        exists || println(io, join(csv_escape.(String.(names)), ','))
        println(io, join((csv_escape(getproperty(row, n)) for n in names), ','))
    end
end

function write_rows(path, rows)
    isempty(rows) && error("no rows for $path")
    open(path, "w") do io
        names = propertynames(first(rows))
        println(io, join(csv_escape.(String.(names)), ','))
        for row in rows
            println(io, join((csv_escape(getproperty(row, n)) for n in names), ','))
        end
    end
end

rho_partial(rho, k) = rho / (1 + k * rho)

function equicorr_matrix(d, rho)
    S = fill(Float64(rho), d, d)
    @inbounds for i in 1:d
        S[i, i] = 1.0
    end
    return S
end

function equicorr_dvine(d, rho)
    levels = [Tuple(GaussianCopula(2, rho_partial(rho, t - 1)) for _ in 1:(d - t)) for t in 1:(d - 1)]
    return DVineCopula(collect(1:d), levels)
end

function equicorr_reference(d, rho, event)
    if event == "cdf"
        lo, hi = zeros(d), fill(1 - 1 / d, d)
    elseif event == "interior"
        lo, hi = fill(1 / (2d), d), fill(1 - 1 / (2d), d)
    else
        error("unknown event $event")
    end
    a = StatsFuns.norminvcdf(hi[1])
    b = StatsFuns.norminvcdf(lo[1])
    sr, so = sqrt(rho), sqrt(1 - rho)
    invsqrt2pi = inv(sqrt(2pi))
    f(z) = begin
        upper = StatsFuns.normcdf((a - sr * z) / so)
        prob = if event == "cdf"
            upper
        else
            lower = StatsFuns.normcdf((b - sr * z) / so)
            max(0.0, upper - lower)
        end
        prob <= 0 && return 0.0
        return exp(-0.5 * z * z + d * log(prob)) * invsqrt2pi
    end
    val, err = Copulas.QuadGK.quadgk(f, -Inf, Inf; rtol=1e-10, atol=1e-12)
    return lo, hi, val, err
end

function rectangle_indicator(U, lower, upper)
    count = 0
    @inbounds for col in axes(U, 2)
        inside = true
        for row in axes(U, 1)
            if !(lower[row] < U[row, col] <= upper[row])
                inside = false
                break
            end
        end
        count += inside
    end
    return count / size(U, 2)
end

function design_points(rule::Symbol, d::Int, N::Int, seed::Int)
    rng = MersenneTwister(seed)
    if rule === :sobol_owen
        return VineCopulas._reduced_owen_points(d, N, UInt64(seed))
    elseif rule === :randomized_halton
        return QuasiMonteCarlo.sample(N, d, QuasiMonteCarlo.RandomizedHaltonSample(rng=rng))
    else
        error("unknown rule $rule")
    end
end

function reduced_from_design!(runner::PublicationReducedRunner, Zred)
    reset_counts!(runner)
    total = 0.0
    @inbounds for col in axes(Zred, 2)
        G, _ = VineCopulas._reduced_probability_integrand!(
            runner.workspace, runner.plan, runner.lower, runner.upper,
            view(Zred, :, col), runner.order, runner.slots, runner.conditioning_values,
        )
        total += G
    end
    return (estimate=total / size(Zred, 2), hfunc_calls=runner.workspace.counts.hfunc,
            hinv_calls=runner.workspace.counts.hinv, pair_cdf_calls=runner.workspace.counts.pair_cdf)
end

function run_cell!(path, case_id, vine, lower, upper, reference, rho, event, d, N, reps, method, rule)
    done = Set((r["case_id"], r["N"], r["rep"], r["method"], r["rule"]) for r in read_csv(path))
    runner = prepare_reduced_runner((vine=vine, lower=lower, upper=upper))
    s = d - 2
    for rep in 1:reps
        key = (case_id, string(N), string(rep), method, string(rule))
        key in done && continue
        seed = BASE_SEED + 1_000_000 * d + 10_000 * round(Int, 100rho) + 100 * N + rep
        gt = @timed design_points(rule, d, N, seed)
        Z = gt.value
        if method == "current_indicator"
            et = @timed rectangle_indicator(inverse_rosenblatt(vine, Z), lower, upper)
            estimate = et.value
            hfunc = NaN; hinv = NaN; paircdf = NaN
        else
            et = @timed reduced_from_design!(runner, view(Z, 1:s, :))
            estimate = et.value.estimate
            hfunc = et.value.hfunc_calls; hinv = et.value.hinv_calls; paircdf = et.value.pair_cdf_calls
        end
        append_row(path, (;
            case_id, d, conditioning_dimension=s, rho, event, rule=string(rule),
            N, rep, seed, method, estimate, reference, abs_error=abs(estimate - reference),
            generation_time_s=gt.time, evaluation_time_s=et.time,
            total_runtime_s=gt.time + et.time, generation_bytes=gt.bytes,
            evaluation_bytes=et.bytes, generation_allocations=Base.gc_alloc_count(gt.gcstats),
            evaluation_allocations=Base.gc_alloc_count(et.gcstats),
            hfunc_calls=hfunc, hinv_calls=hinv, pair_cdf_calls=paircdf,
        ))
    end
end

function projected_seconds(path, case_id, method, rule, N, reps)
    rows = filter(read_csv(path)) do r
        r["case_id"] == case_id && r["method"] == method &&
            r["rule"] == string(rule) && parse(Int, r["N"]) < N
    end
    isempty(rows) && return NaN
    byN = Dict{Int,Vector{Float64}}()
    for r in rows
        n = parse(Int, r["N"])
        push!(get!(byN, n, Float64[]), parse(Float64, r["total_runtime_s"]) / n)
    end
    n0 = maximum(keys(byN))
    return median(byN[n0]) * N * reps
end

function summarize(path)
    rows = read_csv(path)
    groups = Dict{Tuple{String,String,Int,String},Vector{Dict{String,String}}}()
    for r in rows
        push!(get!(groups, (r["case_id"], r["rule"], parse(Int, r["N"]), r["method"]), Dict{String,String}[]), r)
    end
    out = NamedTuple[]
    for key in sort(collect(keys(groups)))
        rr = groups[key]
        errs = [parse(Float64, r["abs_error"]) for r in rr]
        runt = [parse(Float64, r["total_runtime_s"]) for r in rr]
        push!(out, (;
            case_id=key[1], rule=key[2], N=key[3], method=key[4], n=length(rr),
            d=parse(Int, first(rr)["d"]), rho=parse(Float64, first(rr)["rho"]),
            event=first(rr)["event"], rmse=sqrt(mean(abs2, errs)),
            median_abs_error=median(errs), median_total_runtime_s=median(runt),
            runtime_per_point_s=median(runt) / key[3],
        ))
    end
    return out
end

function run_references!()
    existing = Set(r["case_id"] for r in read_csv(REF_OUT))
    for d in (20, 50, 100, 300, 500), rho in (0.3, 0.8), event in ("cdf", "interior")
        case_id = "equicorr_$(event)_rho$(replace(string(rho), "."=>"_"))_d$d"
        case_id in existing && continue
        lo, hi, val, err = equicorr_reference(d, rho, event)
        append_row(REF_OUT, (;
            case_id, d, rho, event, lower=first(lo), upper=first(hi),
            reference=val, quadrature_error=err, reference_method="1D equicorrelated Gaussian mixture quadrature",
        ))
    end
end

function run_density_parity!()
    isfile(PARITY_OUT) && return
    rng = MersenneTwister(BASE_SEED)
    for d in (5, 10), rho in (0.3, 0.8), i in 1:10
        vine = equicorr_dvine(d, rho)
        direct = GaussianCopula(equicorr_matrix(d, rho))
        u = clamp.(rand(rng, d), 0.05, 0.95)
        lv, lg = logpdf(vine, u), logpdf(direct, u)
        append_row(PARITY_OUT, (;
            d, rho, point=i, logpdf_dvine=lv, logpdf_direct=lg,
            abs_difference=abs(lv - lg),
        ))
    end
end

function ref_lookup()
    rows = read_csv(REF_OUT)
    Dict(r["case_id"] => parse(Float64, r["reference"]) for r in rows)
end

run_references!()
run_density_parity!()
refs = ref_lookup()

if !FINALIZE_ONLY
    for d in (20, 50, 100), rho in (0.3, 0.8), event in ("cdf", "interior")
        case_id = "equicorr_$(event)_rho$(replace(string(rho), "."=>"_"))_d$d"
        lo, hi, ref, _ = equicorr_reference(d, rho, event)
        vine = equicorr_dvine(d, rho)
        for N in (256, 1024, 4096, 16384, 65536)
            PILOT_ONLY && N > 256 && continue
            for method in ("current_indicator", "reduced_conditional")
                proj = projected_seconds(RAW_OUT, case_id, method, :sobol_owen, N, PRIMARY_REPS)
                if isfinite(proj) && proj > MAX_CELL_SECONDS
                    append_row(SKIP_OUT, (;
                        case_id, d, rho, event, rule="sobol_owen", N, method,
                        reps=PRIMARY_REPS, projected_seconds=proj,
                        threshold_seconds=MAX_CELL_SECONDS, reason="projected primary cell runtime exceeds threshold",
                    ))
                    continue
                end
                run_cell!(RAW_OUT, case_id, vine, lo, hi, refs[case_id], rho, event, d, N, PRIMARY_REPS, method, :sobol_owen)
            end
        end
    end

    if RUN_STRESS
    for d in (300, 500), rho in (0.3, 0.8), event in ("cdf", "interior")
        case_id = "equicorr_$(event)_rho$(replace(string(rho), "."=>"_"))_d$d"
        lo, hi, ref, _ = equicorr_reference(d, rho, event)
        vine = equicorr_dvine(d, rho)
        for N in (256, 1024, 4096)
            PILOT_ONLY && N > 256 && continue
            for method in ("current_indicator", "reduced_conditional")
                if method == "current_indicator" && SKIP_STRESS_INDICATOR
                    append_row(SKIP_OUT, (;
                        case_id, d, rho, event, rule="sobol_owen", N, method,
                        reps=STRESS_REPS, projected_seconds=NaN,
                        threshold_seconds=MAX_CELL_SECONDS,
                        reason="stress-test full-dimensional indicator skipped after interrupted d=300 N=256 pilot exceeded two minutes inside full inverse Rosenblatt",
                    ))
                    continue
                end
                proj = projected_seconds(STRESS_OUT, case_id, method, :sobol_owen, N, STRESS_REPS)
                if isfinite(proj) && proj > MAX_CELL_SECONDS
                    append_row(SKIP_OUT, (;
                        case_id, d, rho, event, rule="sobol_owen", N, method,
                        reps=STRESS_REPS, projected_seconds=proj,
                        threshold_seconds=MAX_CELL_SECONDS, reason="projected stress cell runtime exceeds threshold",
                    ))
                    continue
                end
                run_cell!(STRESS_OUT, case_id, vine, lo, hi, refs[case_id], rho, event, d, N, STRESS_REPS, method, :sobol_owen)
            end
        end
    end
    end

    if RUN_HALTON_SOBOL
    for d in (50, 100, 300, 500), rho in (0.3, 0.8), event in ("cdf",)
        case_id = "equicorr_$(event)_rho$(replace(string(rho), "."=>"_"))_d$d"
        lo, hi, ref, _ = equicorr_reference(d, rho, event)
        vine = equicorr_dvine(d, rho)
        for N in (256, 1024, 4096), rule in (:sobol_owen, :randomized_halton)
            PILOT_ONLY && N > 256 && continue
            run_cell!(HALTON_SOBOL_OUT, case_id, vine, lo, hi, refs[case_id], rho, event, d, N, 5, "reduced_conditional", rule)
        end
    end
    end
end

summaries = vcat(summarize(RAW_OUT), summarize(STRESS_OUT), summarize(HALTON_SOBOL_OUT))
write_rows(SUMMARY_OUT, summaries)

open(joinpath(OUTDIR, "04_highdim_scaling.md"), "w") do io
    parity = read_csv(PARITY_OUT)
    maxdiff = maximum(parse(Float64, r["abs_difference"]) for r in parity)
    println(io, "# Phase 4 adaptive high-dimensional scaling")
    println(io)
    if FINALIZE_ONLY || PILOT_ONLY
        println(io, "Status: PARTIAL adaptive pilot checkpoint.")
    else
        println(io, "Status: COMPLETE for the adaptive pilot design.")
    end
    println(io)
    println(io, "## Mathematical setup")
    println(io)
    println(io, "An equicorrelated Gaussian copula with common correlation rho is represented as a Gaussian D-vine whose tree-t pair correlations are partial correlations `rho / (1 + (t-1)rho)`. This was validated numerically by comparing the D-vine log density against a direct Gaussian copula for d=5 and d=10.")
    println(io)
    println(io, "- Maximum low-dimensional log-density discrepancy: $maxdiff")
    println(io)
    println(io, "Reference probabilities use the independent one-dimensional Gaussian-mixture representation, not the reduced estimator.")
    println(io)
    println(io, "## Adaptive execution")
    println(io)
    println(io, "- Primary dimensions: d=20,50,100.")
    println(io, "- Stress dimensions: d=300,500.")
    println(io, "- Primary rule: Sobol/Owen.")
    println(io, "- Stress tests and high-dimensional Halton/Sobol pilots are checkpointed separately.")
    println(io, "- Stress section enabled in this run: $RUN_STRESS.")
    println(io, "- Halton/Sobol high-dimensional pilot enabled in this run: $RUN_HALTON_SOBOL.")
    println(io, "- Projection threshold per method cell: $MAX_CELL_SECONDS seconds.")
    println(io, "- Full-dimensional stress indicator skipped by default: $SKIP_STRESS_INDICATOR.")
    println(io, "  An attempted d=300, N=256, one-rep pilot was manually interrupted after more than two minutes inside the full inverse Rosenblatt path, so the stress section records Reduced results and explicit skipped indicator cells instead of fabricating a matched campaign.")
    println(io)
    println(io, "Removing two dimensions is only a small relative dimension reduction at d=300 or d=500; any advantage observed there must come from indicator removal, density-Jacobian cancellation, and integrand regularity, not from dimension reduction alone.")
    println(io)
    println(io, "## Outputs")
    println(io)
    println(io, "- `highdim_reference_values.csv`")
    println(io, "- `highdim_density_parity.csv`")
    println(io, "- `highdim_raw.csv`")
    println(io, "- `highdim_stress_raw.csv`")
    println(io, "- `highdim_halton_sobol.csv`")
    println(io, "- `highdim_summary.csv`")
    println(io, "- `highdim_skipped_cells.csv` if any cells were capped by projected runtime")
end

println("Wrote Phase 4 outputs to $OUTDIR")
