using Copulas
using Distributions
using QuasiMonteCarlo
using Random
using Statistics
using VineCopulas

include(joinpath(@__DIR__, "..", "..", "vendor", "VineCopulas", "benchmarks", "publication", "publication_common.jl"))

const OUTDIR = joinpath(@__DIR__, "rerun")
mkpath(OUTDIR)
const FROZEN = joinpath(@__DIR__, "..", "..", "internal_benchmark", "frozen_results")
const RAW_OUT = joinpath(OUTDIR, "integration_rule_raw.csv")
const SUMMARY_OUT = joinpath(OUTDIR, "integration_rule_summary.csv")
const SKIP_OUT = joinpath(OUTDIR, "integration_rule_skipped_cells.csv")
const POWERS = parse.(Int, split(get(ENV, "POWERS", "8,10,12,14,16"), ','))
const REPS = parse(Int, get(ENV, "REPS", "30"))
const BASE_SEED = parse(Int, get(ENV, "BASE_SEED", "424242"))
const MAX_PROJECTED_CELL_SECONDS = parse(Float64, get(ENV, "MAX_PROJECTED_CELL_SECONDS", "180"))

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

csv_escape(x) = begin
    s = string(x)
    (occursin(',', s) || occursin('"', s) || occursin('\n', s)) ?
        "\"" * replace(s, "\"" => "\"\"") * "\"" : s
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

function read_csv_dicts(path)
    lines = readlines(path)
    header = csv_parse_line(first(lines))
    rows = Vector{Dict{String,String}}()
    for line in Iterators.drop(lines, 1)
        isempty(strip(line)) && continue
        push!(rows, Dict(zip(header, csv_parse_line(line))))
    end
    return rows
end

function references()
    rows = read_csv_dicts(joinpath(FROZEN, "references.csv"))
    return Dict(r["case_id"] => parse(Float64, r["estimate"]) for r in rows)
end

function completed_keys(path)
    isfile(path) || return Set{Tuple{String,String,Int,Int,String}}()
    rows = read_csv_dicts(path)
    return Set((r["case_id"], r["rule"], parse(Int, r["N"]), parse(Int, r["rep"]), r["method"]) for r in rows)
end

function existing_rows(path)
    isfile(path) || return Dict{String,String}[]
    return read_csv_dicts(path)
end

function selected_cases()
    cfg = merge(publication_config(), (; dims=[5, 10, 20], case_set="core"))
    all = publication_cases(cfg)
    ids = Set(["gaussian_moderate_d5", "gaussian_strong_d20", "mixed_d10", "genuine_rvine_d5"])
    return filter(c -> c.id in ids, all)
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
    if rule === :iid_mc
        return rand(rng, d, N)
    elseif rule === :randomized_halton
        return QuasiMonteCarlo.sample(N, d, QuasiMonteCarlo.RandomizedHaltonSample(rng=rng))
    elseif rule === :sobol_owen
        return VineCopulas._reduced_owen_points(d, N, UInt64(seed))
    elseif rule === :shifted_lattice
        return QuasiMonteCarlo.sample(N, d, QuasiMonteCarlo.LatticeRuleSample(R=QuasiMonteCarlo.Shift(rng=rng)))
    else
        error("unknown rule $rule")
    end
end

function reduced_from_design!(runner::PublicationReducedRunner, Zred)
    reset_counts!(runner)
    total = 0.0
    @inbounds for col in axes(Zred, 2)
        G, _ = VineCopulas._reduced_probability_integrand!(
            runner.workspace,
            runner.plan,
            runner.lower,
            runner.upper,
            view(Zred, :, col),
            runner.order,
            runner.slots,
            runner.conditioning_values,
        )
        total += G
    end
    return (;
        estimate=total / size(Zred, 2),
        hfunc_calls=runner.workspace.counts.hfunc,
        hinv_calls=runner.workspace.counts.hinv,
        pair_cdf_calls=runner.workspace.counts.pair_cdf,
    )
end

function summarize_rows(rows)
    groups = Dict{Tuple{String,String,Int,String},Vector{Dict{String,String}}}()
    for r in rows
        push!(get!(groups, (r["case_id"], r["rule"], parse(Int, r["N"]), r["method"]), Dict{String,String}[]), r)
    end
    base = NamedTuple[]
    for key in sort(collect(keys(groups)))
        rr = groups[key]
        errs = [parse(Float64, r["abs_error"]) for r in rr]
        runt = [parse(Float64, r["total_runtime_s"]) for r in rr]
        evalt = [parse(Float64, r["evaluation_time_s"]) for r in rr]
        bias = mean(parse(Float64, r["estimate"]) - parse(Float64, r["reference"]) for r in rr)
        push!(base, (;
            case_id=key[1], rule=key[2], N=key[3], method=key[4], n=length(rr),
            bias=bias, rmse=sqrt(mean(abs2, errs)), median_abs_error=median(errs),
            median_total_runtime_s=median(runt), median_evaluation_time_s=median(evalt),
            p90_abs_error=quantile(errs, 0.9),
        ))
    end
    bycell = Dict((r.case_id, r.rule, r.N) => Dict{String,Any}() for r in base)
    for r in base
        bycell[(r.case_id, r.rule, r.N)][r.method] = r
    end
    out = NamedTuple[]
    for r in base
        cell = bycell[(r.case_id, r.rule, r.N)]
        complete_pair = haskey(cell, "current_indicator") &&
            haskey(cell, "reduced_conditional") &&
            cell["current_indicator"].n == REPS &&
            cell["reduced_conditional"].n == REPS
        if complete_pair
            c = cell["current_indicator"]
            red = cell["reduced_conditional"]
            push!(out, merge(r, (;
                complete_paired_comparison=true,
                rmse_ratio_reduced_current=red.rmse / c.rmse,
                runtime_ratio_reduced_current=red.median_total_runtime_s / c.median_total_runtime_s,
            )))
        else
            push!(out, merge(r, (;
                complete_paired_comparison=false,
                rmse_ratio_reduced_current=NaN,
                runtime_ratio_reduced_current=NaN,
            )))
        end
    end
    return out
end

function target_powers(case_id)
    powers = copy(POWERS)
    case_id == "gaussian_strong_d20" && filter!(<=(14), powers)
    return powers
end

function projected_cell_seconds(rows, case_id, rule, method, N, remaining_reps)
    candidates = filter(rows) do r
        r["case_id"] == case_id &&
            r["rule"] == string(rule) &&
            r["method"] == method &&
            parse(Int, r["N"]) < N
    end
    isempty(candidates) && return NaN
    byN = Dict{Int,Vector{Float64}}()
    for r in candidates
        n = parse(Int, r["N"])
        push!(get!(byN, n, Float64[]), parse(Float64, r["total_runtime_s"]) / n)
    end
    n0 = maximum(keys(byN))
    per_point = median(byN[n0])
    return per_point * N * remaining_reps
end

function append_skip(row)
    append_row(SKIP_OUT, row)
end

rules = [:iid_mc, :randomized_halton, :sobol_owen, :shifted_lattice]
refs = references()
done = completed_keys(RAW_OUT)
rows_seen = existing_rows(RAW_OUT)

for case in selected_cases()
    runner = prepare_reduced_runner(case)
    d = length(case.vine)
    s = d - 2
    # Warm up both methods and all rules outside timing.
    Zwarm = design_points(:sobol_owen, d, 16, BASE_SEED)
    rectangle_indicator(inverse_rosenblatt(case.vine, Zwarm), case.lower, case.upper)
    reduced_from_design!(runner, view(Zwarm, 1:s, :))

    for rule in rules, pow in target_powers(case.id)
        N = 1 << pow
        for method in ("current_indicator", "reduced_conditional")
            remaining = count(rep -> !((case.id, string(rule), N, rep, method) in done), 1:REPS)
            if remaining > 0
                projected = projected_cell_seconds(rows_seen, case.id, rule, method, N, remaining)
                if isfinite(projected) && projected > MAX_PROJECTED_CELL_SECONDS
                    append_skip((;
                        case_id=case.id, rule=string(rule), N, method,
                        remaining_reps=remaining, projected_seconds=projected,
                        threshold_seconds=MAX_PROJECTED_CELL_SECONDS,
                        reason="projected cell runtime exceeds threshold",
                    ))
                    continue
                end
            end
            for rep in 1:REPS
            key = (case.id, string(rule), N, rep, method)
            key in done && continue
            seed = BASE_SEED + 100000 * findfirst(==(rule), rules) + 1000 * pow + rep
            gt = @timed design_points(rule, d, N, seed)
            Z = gt.value
            if method == "current_indicator"
                et = @timed begin
                    U = inverse_rosenblatt(case.vine, Z)
                    rectangle_indicator(U, case.lower, case.upper)
                end
                estimate = et.value
                hfunc = NaN
                hinv = NaN
                paircdf = NaN
            else
                et = @timed reduced_from_design!(runner, view(Z, 1:s, :))
                estimate = et.value.estimate
                hfunc = et.value.hfunc_calls
                hinv = et.value.hinv_calls
                paircdf = et.value.pair_cdf_calls
            end
            ref = refs[case.id]
            append_row(RAW_OUT, (;
                case_id=case.id, design=case.design, model=case.model, point=case.point,
                d, conditioning_dimension=s, rule=string(rule), N, rep, seed,
                method, estimate, reference=ref, abs_error=abs(estimate - ref),
                generation_time_s=gt.time, evaluation_time_s=et.time,
                total_runtime_s=gt.time + et.time,
                generation_bytes=gt.bytes, evaluation_bytes=et.bytes,
                generation_allocations=Base.gc_alloc_count(gt.gcstats),
                evaluation_allocations=Base.gc_alloc_count(et.gcstats),
                hfunc_calls=hfunc, hinv_calls=hinv, pair_cdf_calls=paircdf,
            ))
            push!(done, key)
            push!(rows_seen, Dict(
                "case_id" => case.id, "rule" => string(rule), "N" => string(N),
                "method" => method, "total_runtime_s" => string(gt.time + et.time),
            ))
            end
        end
    end
end

summary = summarize_rows(read_csv_dicts(RAW_OUT))
write_rows(SUMMARY_OUT, summary)
complete_pairs = filter(r -> r.method == "reduced_conditional" && r.complete_paired_comparison, summary)
full_budget_text = join((string(1 << p) for p in POWERS), ", ")

open(joinpath(OUTDIR, "03_integration_rule_sensitivity.md"), "w") do io
    println(io, "# Phase 3 integration-rule sensitivity")
    println(io)
    println(io, "Status: COMPLETE under the revised adaptive design.")
    println(io)
    println(io, "This experiment compares the same full-dimensional indicator representation and the same reduced representation under several point designs. It does not change the mathematical reduction.")
    println(io)
    println(io, "## Design")
    println(io)
    println(io, "- Cases: gaussian_moderate_d5, gaussian_strong_d20, mixed_d10, genuine_rvine_d5")
    println(io, "- Budgets for gaussian_moderate_d5, mixed_d10, genuine_rvine_d5: $full_budget_text")
    println(io, "- Budgets for gaussian_strong_d20 paired comparison: 256, 1024, 4096, 16384")
    println(io, "- Replicates per cell: $REPS")
    println(io, "- Rules: iid MC, randomized Halton, Sobol with Owen nested uniform scrambling, shifted rank-1 lattice")
    println(io, "- Coupling: each replicate generates a d-dimensional design once; Indicator uses all d coordinates and Reduced uses the first d-2 coordinates.")
    println(io)
    println(io, "## Lattice and Halton qualification")
    println(io)
    println(io, "`RandomizedHaltonSample` is documented in the installed QuasiMonteCarlo.jl 0.3.11 and cites Owen's randomized Halton algorithm. `LatticeRuleSample(R=Shift(...))` uses the package's LatticeRules.jl-backed rank-1 lattice with Cranley-Patterson shift. No scrambled-net theorem is claimed for these rules.")
    println(io)
    println(io, "## Outputs")
    println(io)
    println(io, "- Raw replicate data: `integration_rule_raw.csv`")
    println(io, "- Aggregated summary: `integration_rule_summary.csv`")
    println(io, "- Skipped/projection log if any: `integration_rule_skipped_cells.csv`")
    println(io)
    println(io, "Raw rows retained: $(length(read_csv_dicts(RAW_OUT))).")
    if isfile(SKIP_OUT)
        skipped = read_csv_dicts(SKIP_OUT)
        println(io, "Projected-cost skips: $(length(skipped)).")
    else
        println(io, "Projected-cost skips: 0.")
    end
    println(io)
    println(io, "## Paired complete-cell summary")
    println(io)
    println(io, "- Complete paired reduced/current cells: $(length(complete_pairs))")
    if !isempty(complete_pairs)
        ratios = sort([r.rmse_ratio_reduced_current for r in complete_pairs])
        runtimes = sort([r.runtime_ratio_reduced_current for r in complete_pairs])
        println(io, "- Median RMSE ratio reduced/current: $(median(ratios))")
        println(io, "- Median runtime ratio reduced/current: $(median(runtimes))")
        println(io, "- Complete cells with reduced RMSE lower than current: $(count(<(1.0), ratios))/$(length(ratios))")
    end
    println(io)
    println(io, "The d=20 full-indicator path is substantially more expensive because it evaluates the full inverse Rosenblatt transform before applying the rectangle indicator. The d=20 campaign is intentionally capped at N=16384 for paired comparisons; any previously generated N=65536 rows are retained as secondary incomplete data and are not used for reduced/current ratio columns.")
    println(io)
    println(io, "Figures are generated by `phase3_make_figures.jl` from the retained CSV files.")
end

println("Wrote Phase 3 integration-rule sensitivity outputs to $OUTDIR")
