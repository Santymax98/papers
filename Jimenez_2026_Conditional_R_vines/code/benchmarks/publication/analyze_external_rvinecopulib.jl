include("publication_common.jl")

const CFG = publication_config()
const EXT_DIR = joinpath(CFG.out, "external_rvinecopulib")
const EXT_RAW_PATH = joinpath(EXT_DIR, "raw_results.csv")
const INT_SUMMARY_PATH = joinpath(CFG.out, "summary_by_case.csv")
const REF_PATH = joinpath(CFG.out, "references.csv")
const CASES_PATH = joinpath(EXT_DIR, "cases.csv")
const AUDIT_PATH = joinpath(EXT_DIR, "model_parity_audit.csv")

ext_raw = read_table(EXT_RAW_PATH)
int_summary_raw = read_table(INT_SUMMARY_PATH)
refs = read_table(REF_PATH)
cases_raw = read_table(CASES_PATH)
audit_raw = read_table(AUDIT_PATH)

isempty(ext_raw) && error("external rvinecopulib raw_results.csv is empty")
isempty(int_summary_raw) && error("summary_by_case.csv is missing; run the internal publication analysis first")
isempty(refs) && error("references.csv is missing")
isempty(cases_raw) && error("external cases.csv is missing")
isempty(audit_raw) && error("model parity audit is missing")
all(_parse_bool(r, "passed") for r in audit_raw) || error("model parity audit contains failures")

ref_by_case = Dict(r["case_id"] => r for r in refs)
case_by_id = Dict(r["case_id"] => r for r in cases_raw)
external_case_ids = Set(keys(case_by_id))

rmse(xs) = isempty(xs) ? NaN : sqrt(mean(abs2, xs))

function bootstrap_ci(values::Vector{Float64}, stat, B::Int, rng::AbstractRNG; alpha=0.05)
    n = length(values)
    n == 0 && return (NaN, NaN)
    n == 1 && return (stat(values), stat(values))
    boot = Vector{Float64}(undef, B)
    idx = Vector{Int}(undef, n)
    @inbounds for b in 1:B
        for j in 1:n
            idx[j] = rand(rng, 1:n)
        end
        boot[b] = stat(values[idx])
    end
    return quantile(boot, alpha / 2), quantile(boot, 1 - alpha / 2)
end

# -----------------------------------------------------------------------------
# External cell summaries, using the same reference-resolution and bootstrap
# policy as the internal publication analysis.
# -----------------------------------------------------------------------------

groups = Dict{Tuple{String,Int},Vector{Dict{String,String}}}()
for row in ext_raw
    key = (row["case_id"], _parse_int(row, "N"))
    push!(get!(groups, key, Dict{String,String}[]), row)
end

external_summary = NamedTuple[]
for (group_idx, key) in enumerate(sort!(collect(keys(groups)); by=string))
    rows = groups[key]
    errors = Float64[_parse_float(r, "abs_error") for r in rows]
    runtimes = Float64[_parse_float(r, "runtime_s") for r in rows]
    estimates = Float64[_parse_float(r, "estimate") for r in rows]
    firstrow = first(rows)
    ref = ref_by_case[key[1]]
    ref_unc = _parse_float(ref, "uncertainty")
    ref_split_z = _parse_float(ref, "split_half_z")
    value_rmse = rmse(errors)
    ratio = value_rmse > 0 ? ref_unc / value_rmse : Inf
    split_stable = isfinite(ref_split_z) && ref_split_z <= CFG.reference_split_z_max
    resolved = split_stable && ratio <= CFG.max_reference_ratio

    rng1 = MersenneTwister(CFG.bootstrap_seed + 1_000_000 + 10 * group_idx)
    rng2 = MersenneTwister(CFG.bootstrap_seed + 1_000_001 + 10 * group_idx)
    rmse_lo, rmse_hi = bootstrap_ci(errors, rmse, CFG.bootstrap_b, rng1)
    rt_lo, rt_hi = bootstrap_ci(runtimes, median, CFG.bootstrap_b, rng2)

    push!(external_summary, (;
        case_id=key[1],
        design=firstrow["design"],
        model=firstrow["model"],
        point=firstrow["point"],
        d=_parse_int(firstrow, "d"),
        N=key[2],
        method="rvinecopulib_pvinecop",
        n=length(rows),
        reference=_parse_float(firstrow, "reference"),
        reference_uncertainty=ref_unc,
        reference_target_met=_parse_bool(ref, "target_met"),
        reference_split_half_z=ref_split_z,
        reference_split_stable=split_stable,
        reference_uncertainty_to_rmse=ratio,
        reference_resolved=resolved,
        mean_estimate=mean(estimates),
        estimate_sd=length(estimates) > 1 ? std(estimates) : 0.0,
        median_abs_error=median(errors),
        mean_abs_error=mean(errors),
        rmse=value_rmse,
        rmse_ci_low=rmse_lo,
        rmse_ci_high=rmse_hi,
        p90_abs_error=quantile(errors, 0.90),
        max_abs_error=maximum(errors),
        median_runtime_s=median(runtimes),
        runtime_ci_low_s=rt_lo,
        runtime_ci_high_s=rt_hi,
    ))
end
sort!(external_summary; by=r -> (r.d, r.case_id, r.N))

external_reps_expected = _env_int("EXTERNAL_REPS", CFG.reps)
expected_N = Set(1 << p for p in CFG.powers)
for case_id in external_case_ids
    rows = filter(r -> r.case_id == case_id, external_summary)
    Set(getfield.(rows, :N)) == expected_N || error("incomplete external N grid for $case_id")
    all(r -> r.n == external_reps_expected, rows) || error("incomplete external replicates for $case_id")
end

write_rows(joinpath(EXT_DIR, "summary.csv"), external_summary)

# -----------------------------------------------------------------------------
# Common three-method summary for CDF cases only.
# -----------------------------------------------------------------------------

combined = NamedTuple[]
for r in int_summary_raw
    r["case_id"] in external_case_ids || continue
    r["method"] in ("current_indicator", "reduced_conditional") || continue
    push!(combined, (;
        case_id=r["case_id"], design=r["design"], model=r["model"], point=r["point"],
        d=_parse_int(r, "d"), N=_parse_int(r, "N"), method=r["method"], n=_parse_int(r, "n"),
        reference=_parse_float(r, "reference"), reference_uncertainty=_parse_float(r, "reference_uncertainty"),
        reference_target_met=_parse_bool(r, "reference_target_met"),
        reference_uncertainty_to_rmse=_parse_float(r, "reference_uncertainty_to_rmse"),
        reference_resolved=_parse_bool(r, "reference_resolved"),
        rmse=_parse_float(r, "rmse"), rmse_ci_low=_parse_float(r, "rmse_ci_low"),
        rmse_ci_high=_parse_float(r, "rmse_ci_high"), median_abs_error=_parse_float(r, "median_abs_error"),
        median_runtime_s=_parse_float(r, "median_runtime_s"), runtime_ci_low_s=_parse_float(r, "runtime_ci_low_s"),
        runtime_ci_high_s=_parse_float(r, "runtime_ci_high_s"),
    ))
end
for r in external_summary
    push!(combined, (;
        case_id=r.case_id, design=r.design, model=r.model, point=r.point,
        d=r.d, N=r.N, method=r.method, n=r.n,
        reference=r.reference, reference_uncertainty=r.reference_uncertainty,
        reference_target_met=r.reference_target_met,
        reference_uncertainty_to_rmse=r.reference_uncertainty_to_rmse,
        reference_resolved=r.reference_resolved,
        rmse=r.rmse, rmse_ci_low=r.rmse_ci_low, rmse_ci_high=r.rmse_ci_high,
        median_abs_error=r.median_abs_error,
        median_runtime_s=r.median_runtime_s, runtime_ci_low_s=r.runtime_ci_low_s,
        runtime_ci_high_s=r.runtime_ci_high_s,
    ))
end
sort!(combined; by=r -> (r.d, r.case_id, r.N, r.method))
write_rows(joinpath(EXT_DIR, "three_method_summary.csv"), combined)

lookup = Dict((r.case_id, r.N, r.method) => r for r in combined)

# Equal-N comparisons are cell-level only. The R and Julia randomized designs
# are independent; we intentionally do NOT claim replicate-level pairing across
# implementations.
equal_budget = NamedTuple[]
for ext in external_summary
    rk = (ext.case_id, ext.N, "reduced_conditional")
    ck = (ext.case_id, ext.N, "current_indicator")
    haskey(lookup, rk) && haskey(lookup, ck) || continue
    red = lookup[rk]
    cur = lookup[ck]
    resolved = ext.reference_resolved && red.reference_resolved && cur.reference_resolved
    push!(equal_budget, (;
        case_id=ext.case_id,
        d=ext.d,
        N=ext.N,
        all_reference_resolved=resolved,
        current_rmse=cur.rmse,
        rvinecopulib_rmse=ext.rmse,
        reduced_rmse=red.rmse,
        reduced_over_rvinecopulib_rmse=ext.rmse > 0 ? red.rmse / ext.rmse : Inf,
        current_over_rvinecopulib_rmse=ext.rmse > 0 ? cur.rmse / ext.rmse : Inf,
        current_runtime_s=cur.median_runtime_s,
        rvinecopulib_runtime_s=ext.median_runtime_s,
        reduced_runtime_s=red.median_runtime_s,
        reduced_over_rvinecopulib_runtime=ext.median_runtime_s > 0 ? red.median_runtime_s / ext.median_runtime_s : Inf,
        reduced_speedup_vs_rvinecopulib=red.median_runtime_s > 0 ? ext.median_runtime_s / red.median_runtime_s : Inf,
        reduced_lower_rmse=red.rmse < ext.rmse,
        reduced_faster=red.median_runtime_s < ext.median_runtime_s,
        reduced_wins_both=(red.rmse < ext.rmse && red.median_runtime_s < ext.median_runtime_s),
        rvinecopulib_lower_rmse_than_current=ext.rmse < cur.rmse,
        rvinecopulib_faster_than_current=ext.median_runtime_s < cur.median_runtime_s,
    ))
end
sort!(equal_budget; by=r -> (r.d, r.case_id, r.N))
write_rows(joinpath(EXT_DIR, "equal_budget.csv"), equal_budget)

# -----------------------------------------------------------------------------
# Matched accuracy: every resolved rvinecopulib target is retained. For each
# target, choose the FASTEST reduced cell whose RMSE is no larger. No
# cherry-picking of the largest observed speedup per case.
# -----------------------------------------------------------------------------

by_case = Dict{String,Vector{typeof(first(combined))}}()
for r in combined
    push!(get!(by_case, r.case_id, typeof(r)[]), r)
end

matched = NamedTuple[]
strictest = NamedTuple[]
for case_id in sort(collect(keys(by_case)))
    rows = by_case[case_id]
    targets = sort(filter(r -> r.method == "rvinecopulib_pvinecop" && r.reference_resolved, rows); by=r -> r.N)
    reduced = filter(r -> r.method == "reduced_conditional" && r.reference_resolved, rows)
    isempty(targets) && continue
    isempty(reduced) && continue

    case_rows = NamedTuple[]
    for target in targets
        candidates = filter(r -> r.rmse <= target.rmse, reduced)
        if isempty(candidates)
            row = (;
                case_id, d=target.d, target_rvinecopulib_N=target.N,
                target_rmse=target.rmse, target_runtime_s=target.median_runtime_s,
                matched=false, reduced_N=0, reduced_rmse=NaN,
                reduced_runtime_s=NaN, speedup=NaN, N_ratio=NaN, rmse_ratio=NaN,
            )
        else
            candidate = candidates[argmin(getfield.(candidates, :median_runtime_s))]
            row = (;
                case_id, d=target.d, target_rvinecopulib_N=target.N,
                target_rmse=target.rmse, target_runtime_s=target.median_runtime_s,
                matched=true, reduced_N=candidate.N, reduced_rmse=candidate.rmse,
                reduced_runtime_s=candidate.median_runtime_s,
                speedup=target.median_runtime_s / candidate.median_runtime_s,
                N_ratio=target.N / candidate.N,
                rmse_ratio=candidate.rmse / target.rmse,
            )
        end
        push!(matched, row)
        push!(case_rows, row)
    end

    strict_target = targets[argmin(getfield.(targets, :rmse))]
    selected = filter(r -> r.target_rvinecopulib_N == strict_target.N, case_rows)
    isempty(selected) || push!(strictest, only(selected))
end
write_rows(joinpath(EXT_DIR, "matched_accuracy.csv"), matched)
!isempty(strictest) && write_rows(joinpath(EXT_DIR, "matched_accuracy_strictest.csv"), strictest)

# Three-method Pareto frontier by case.
pareto = NamedTuple[]
for case_id in sort(collect(keys(by_case)))
    rows = by_case[case_id]
    for method in ("current_indicator", "rvinecopulib_pvinecop", "reduced_conditional")
        selected = sort(filter(r -> r.method == method && r.reference_resolved, rows); by=r -> r.median_runtime_s)
        best = Inf
        for r in selected
            if r.rmse < best
                push!(pareto, (;
                    case_id, d=r.d, method, N=r.N, rmse=r.rmse,
                    median_runtime_s=r.median_runtime_s,
                    rmse_ci_low=r.rmse_ci_low, rmse_ci_high=r.rmse_ci_high,
                    runtime_ci_low_s=r.runtime_ci_low_s, runtime_ci_high_s=r.runtime_ci_high_s,
                ))
                best = r.rmse
            end
        end
    end
end
!isempty(pareto) && write_rows(joinpath(EXT_DIR, "pareto_frontier.csv"), pareto)

figdir = joinpath(EXT_DIR, "figure_data")
mkpath(figdir)
write_rows(joinpath(figdir, "rmse_vs_N_three_methods.csv"), combined)
write_rows(joinpath(figdir, "rmse_vs_runtime_three_methods.csv"), combined)
write_rows(joinpath(figdir, "equal_budget.csv"), equal_budget)
write_rows(joinpath(figdir, "matched_accuracy.csv"), matched)
!isempty(strictest) && write_rows(joinpath(figdir, "matched_accuracy_strictest.csv"), strictest)

tabledir = joinpath(EXT_DIR, "table_data")
mkpath(tabledir)
write_rows(joinpath(tabledir, "equal_budget.csv"), equal_budget)
write_rows(joinpath(tabledir, "matched_accuracy.csv"), matched)
!isempty(strictest) && write_rows(joinpath(tabledir, "matched_accuracy_strictest.csv"), strictest)

valid_equal = filter(r -> r.all_reference_resolved, equal_budget)
valid_match = filter(r -> r.matched, matched)
valid_strict = filter(r -> r.matched, strictest)
rmse_ratios = [r.reduced_over_rvinecopulib_rmse for r in valid_equal if isfinite(r.reduced_over_rvinecopulib_rmse)]
speedups = [r.speedup for r in valid_match if isfinite(r.speedup)]
strict_speedups = [r.speedup for r in valid_strict if isfinite(r.speedup)]

open(joinpath(EXT_DIR, "REPORT.md"), "w") do io
    println(io, "# External rvinecopulib benchmark report")
    println(io)
    println(io, "## Design")
    println(io)
    println(io, "- CDF cases: $(length(external_case_ids))")
    println(io, "- External method: `rvinecopulib::pvinecop`")
    println(io, "- Budgets: $(join(sort(unique(getfield.(external_summary, :N))), ", "))")
    println(io, "- Replicates per external cell: $(minimum(getfield.(external_summary, :n)))")
    println(io, "- Cores: 1")
    println(io, "- General rectangles are excluded; no 2^d inclusion-exclusion baseline is manufactured.")
    println(io, "- R and Julia models are audited by deterministic log-density equality before timing.")
    println(io, "- Cross-language randomized replicates are independent; equal-budget comparisons below are cell-level, not replicate-paired.")
    println(io)
    println(io, "## Model parity")
    println(io)
    println(io, "- Density audit rows passed: $(count(r -> _parse_bool(r, "passed"), audit_raw))/$(length(audit_raw))")
    println(io)
    println(io, "## Reference resolution")
    println(io)
    println(io, "- External aggregate cells: $(length(external_summary))")
    println(io, "- Resolved external cells: $(count(r -> r.reference_resolved, external_summary))/$(length(external_summary))")
    println(io)
    println(io, "## Equal QMC budget: reduced vs rvinecopulib")
    println(io)
    println(io, "- Resolved three-method cells: $(length(valid_equal))/$(length(equal_budget))")
    println(io, "- Reduced lower RMSE: $(count(r -> r.reduced_lower_rmse, valid_equal))/$(length(valid_equal))")
    println(io, "- Reduced lower median runtime: $(count(r -> r.reduced_faster, valid_equal))/$(length(valid_equal))")
    println(io, "- Reduced wins both: $(count(r -> r.reduced_wins_both, valid_equal))/$(length(valid_equal))")
    !isempty(rmse_ratios) && println(io, "- Median reduced/rvinecopulib RMSE ratio: $(median(rmse_ratios))")
    println(io, "- rvinecopulib lower RMSE than current indicator: $(count(r -> r.rvinecopulib_lower_rmse_than_current, valid_equal))/$(length(valid_equal))")
    println(io)
    println(io, "## Matched accuracy: rvinecopulib targets")
    println(io)
    println(io, "Every resolved rvinecopulib target cell is retained; the fastest reduced cell with RMSE no larger is selected.")
    println(io)
    println(io, "- Target cells: $(length(matched))")
    println(io, "- Targets matched by reduced: $(length(valid_match))/$(length(matched))")
    if !isempty(speedups)
        println(io, "- Median speedup: $(median(speedups))")
        println(io, "- 10th percentile speedup: $(quantile(speedups, 0.10))")
        println(io, "- 90th percentile speedup: $(quantile(speedups, 0.90))")
    end
    if !isempty(strict_speedups)
        println(io, "- Median speedup at the strictest rvinecopulib target per case: $(median(strict_speedups))")
    end
    println(io)
    println(io, "Raw external replicates, model-parity audit, bootstrap summaries, matched-accuracy data, Pareto frontiers, and figure/table-ready CSVs are retained under this directory.")
end

println("External rvinecopulib analysis written to ", EXT_DIR)
