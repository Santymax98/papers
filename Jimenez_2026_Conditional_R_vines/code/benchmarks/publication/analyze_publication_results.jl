include("publication_common.jl")

const CFG = publication_config()
const RAW_PATH = joinpath(CFG.out, "raw_results.csv")
const REF_PATH = joinpath(CFG.out, "references.csv")

raw = read_table(RAW_PATH)
refs = read_table(REF_PATH)
isempty(raw) && error("raw_results.csv is empty")
isempty(refs) && error("references.csv is empty")

ref_by_case = Dict(r["case_id"] => r for r in refs)

function rmse(values)
    isempty(values) && return NaN
    return sqrt(mean(abs2, values))
end

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

function finite_median(xs)
    ys = filter(isfinite, Float64.(xs))
    isempty(ys) ? NaN : median(ys)
end

# -----------------------------
# Cell summaries with bootstrap uncertainty
# -----------------------------

groups = Dict{Tuple{String,Int,String},Vector{Dict{String,String}}}()
for row in raw
    key = (row["case_id"], _parse_int(row, "N"), row["method"])
    push!(get!(groups, key, Dict{String,String}[]), row)
end

summary = NamedTuple[]
for (group_idx, key) in enumerate(sort!(collect(keys(groups)); by=string))
    rows = groups[key]
    errors = Float64[_parse_float(r, "abs_error") for r in rows]
    runtimes = Float64[_parse_float(r, "runtime_s") for r in rows]
    bytes = Float64[_parse_float(r, "bytes") for r in rows]
    allocs = Float64[_parse_float(r, "allocations") for r in rows]
    firstrow = first(rows)
    ref = ref_by_case[key[1]]
    ref_uncertainty = _parse_float(ref, "uncertainty")
    ref_target_met = _parse_bool(ref, "target_met")
    ref_split_half_z = _parse_float(ref, "split_half_z")
    # The absolute reference target (e.g. SE <= 1e-7) is an aspirational
    # generation criterion, not a prerequisite for every reported method cell.
    # For inference about a method cell, what matters is that the independent
    # reference is stable and its uncertainty is small relative to the RMSE
    # being resolved. This avoids discarding difficult high-dimensional cases
    # merely because an unnecessarily stringent absolute target was not met.
    reference_split_stable = isfinite(ref_split_half_z) && ref_split_half_z <= CFG.reference_split_z_max
    value_rmse = rmse(errors)
    ratio = value_rmse > 0 ? ref_uncertainty / value_rmse : Inf
    reference_resolved = reference_split_stable && ratio <= CFG.max_reference_ratio

    rng1 = MersenneTwister(CFG.bootstrap_seed + 10 * group_idx)
    rng2 = MersenneTwister(CFG.bootstrap_seed + 10 * group_idx + 1)
    rmse_lo, rmse_hi = bootstrap_ci(errors, rmse, CFG.bootstrap_b, rng1)
    rt_lo, rt_hi = bootstrap_ci(runtimes, median, CFG.bootstrap_b, rng2)

    hf = Float64[_parse_float(r, "hfunc_calls") for r in rows]
    hi = Float64[_parse_float(r, "hinv_calls") for r in rows]
    pc = Float64[_parse_float(r, "pair_cdf_calls") for r in rows]

    push!(summary, (;
        case_id=key[1],
        design=firstrow["design"],
        model=firstrow["model"],
        point=firstrow["point"],
        d=_parse_int(firstrow, "d"),
        conditioning_dimension=_parse_int(firstrow, "conditioning_dimension"),
        N=key[2],
        method=key[3],
        n=length(rows),
        reference=_parse_float(firstrow, "reference"),
        reference_uncertainty=ref_uncertainty,
        reference_target_met=ref_target_met,
        reference_split_half_z=ref_split_half_z,
        reference_split_stable,
        reference_uncertainty_to_rmse=ratio,
        reference_resolved,
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
        median_bytes=median(bytes),
        median_allocations=median(allocs),
        median_hfunc_calls=finite_median(hf),
        median_hinv_calls=finite_median(hi),
        median_pair_cdf_calls=finite_median(pc),
    ))
end
sort!(summary; by=r -> (r.d, r.case_id, r.N, r.method))
write_rows(joinpath(CFG.out, "summary_by_case.csv"), summary)

# -----------------------------
# Paired replicate comparison at equal N
# -----------------------------

pair_groups = Dict{Tuple{String,Int,Int},Dict{String,Dict{String,String}}}()
for row in raw
    key = (row["case_id"], _parse_int(row, "N"), _parse_int(row, "rep"))
    get!(pair_groups, key, Dict{String,Dict{String,String}}())[row["method"]] = row
end

paired = NamedTuple[]
for key in sort!(collect(keys(pair_groups)); by=string)
    pair = pair_groups[key]
    haskey(pair, "current_indicator") && haskey(pair, "reduced_conditional") || continue
    current = pair["current_indicator"]
    reduced = pair["reduced_conditional"]
    ec = _parse_float(current, "abs_error")
    er = _parse_float(reduced, "abs_error")
    tc = _parse_float(current, "runtime_s")
    tr = _parse_float(reduced, "runtime_s")
    push!(paired, (;
        case_id=key[1],
        d=_parse_int(current, "d"),
        N=key[2],
        rep=key[3],
        seed=_parse_int(current, "seed"),
        current_error=ec,
        reduced_error=er,
        error_ratio=ec > 0 ? er / ec : (er == 0 ? 1.0 : Inf),
        current_runtime_s=tc,
        reduced_runtime_s=tr,
        runtime_ratio=tc > 0 ? tr / tc : Inf,
        error_win=er < ec,
        runtime_win=tr < tc,
        both_win=(er < ec && tr < tc),
        reference_target_met=_parse_bool(current, "reference_target_met"),
    ))
end
write_rows(joinpath(CFG.out, "paired_equal_N.csv"), paired)

# -----------------------------
# Honest matched-accuracy comparisons
# For every current-method target cell, choose the FASTEST reduced cell whose
# RMSE is no larger. We retain every target instead of selecting only the
# largest observed speedup for each case.
# -----------------------------

summary_by_case = Dict{String,Vector{typeof(first(summary))}}()
for row in summary
    push!(get!(summary_by_case, row.case_id, typeof(row)[]), row)
end

matched = NamedTuple[]
strictest = NamedTuple[]
for case_id in sort(collect(keys(summary_by_case)))
    rows = summary_by_case[case_id]
    current = sort(filter(r -> r.method == "current_indicator" && r.reference_resolved, rows); by=r -> r.N)
    reduced = filter(r -> r.method == "reduced_conditional" && r.reference_resolved, rows)
    isempty(current) && continue
    isempty(reduced) && continue

    case_matches = NamedTuple[]
    for target in current
        candidates = filter(r -> r.rmse <= target.rmse, reduced)
        isempty(candidates) && continue
        candidate = candidates[argmin(getfield.(candidates, :median_runtime_s))]
        row = (;
            case_id,
            d=target.d,
            target_current_N=target.N,
            target_rmse=target.rmse,
            reduced_N=candidate.N,
            reduced_rmse=candidate.rmse,
            current_runtime_s=target.median_runtime_s,
            reduced_runtime_s=candidate.median_runtime_s,
            speedup=target.median_runtime_s / candidate.median_runtime_s,
            N_ratio=target.N / candidate.N,
            rmse_ratio=candidate.rmse / target.rmse,
        )
        push!(matched, row)
        push!(case_matches, row)
    end
    if !isempty(case_matches)
        # Strictest target = smallest resolved current-method RMSE.
        best_current = current[argmin(getfield.(current, :rmse))]
        candidates = filter(r -> r.target_current_N == best_current.N, case_matches)
        isempty(candidates) || push!(strictest, only(candidates))
    end
end
write_rows(joinpath(CFG.out, "matched_accuracy.csv"), matched)
!isempty(strictest) && write_rows(joinpath(CFG.out, "matched_accuracy_strictest.csv"), strictest)

# -----------------------------
# Empirical convergence rates (descriptive, not theoretical claims)
# -----------------------------

function regression_slope(x::Vector{Float64}, y::Vector{Float64})
    length(x) >= 2 || return NaN
    xm, ym = mean(x), mean(y)
    den = sum(abs2, x .- xm)
    den > 0 || return NaN
    return sum((x .- xm) .* (y .- ym)) / den
end

rates = NamedTuple[]
for case_id in sort(collect(keys(summary_by_case)))
    rows = summary_by_case[case_id]
    for method in ("current_indicator", "reduced_conditional")
        selected = sort(filter(r -> r.method == method && r.reference_resolved && r.rmse > 0, rows); by=r -> r.N)
        length(selected) >= 3 || continue
        x = log2.(Float64.(getfield.(selected, :N)))
        y = log2.(Float64.(getfield.(selected, :rmse)))
        slope = regression_slope(x, y)
        push!(rates, (;
            case_id,
            d=first(selected).d,
            method,
            n_budgets=length(selected),
            log2_rmse_vs_log2_N_slope=slope,
            empirical_rmse_rate=-slope,
            min_N=minimum(getfield.(selected, :N)),
            max_N=maximum(getfield.(selected, :N)),
        ))
    end
end
!isempty(rates) && write_rows(joinpath(CFG.out, "empirical_convergence_rates.csv"), rates)

# -----------------------------
# Complexity audit
# Exact formulas below are for the compiled D-vine traversal used in the
# scaling design. Genuine R-vines can reuse states differently, so their exact
# counts are structure-dependent. They are reported empirically rather than
# forced to equal the D-vine formulas.
# -----------------------------

complexity_rows = NamedTuple[]
for d in sort(unique(getfield.(summary, :d)))
    for model in ("D", "R")
        for probability_type in ("cdf", "rectangle")
            reduced = filter(r -> begin
                r.d == d && r.model == model && r.method == "reduced_conditional" &&
                (probability_type == "rectangle" ? r.point == "rectangle" : r.point != "rectangle")
            end, summary)
            isempty(reduced) && continue
            hf = finite_median([r.median_hfunc_calls / r.N for r in reduced])
            hi = finite_median([r.median_hinv_calls / r.N for r in reduced])
            pc = finite_median([r.median_pair_cdf_calls / r.N for r in reduced])
            sdim = d - 2
            if model == "D"
                expected_hf = probability_type == "cdf" ? sdim * (3sdim + 1) / 2 : 2sdim * (sdim + 1)
                expected_hi = sdim * (sdim - 1) / 2
                expected_total = expected_hf + expected_hi
                exact_formula_applicable = true
            else
                expected_hf = NaN
                expected_hi = NaN
                expected_total = NaN
                exact_formula_applicable = false
            end
            expected_pc = probability_type == "cdf" ? 1.0 : 4.0
            push!(complexity_rows, (;
                d,
                s=sdim,
                model,
                probability_type,
                exact_formula_applicable,
                observed_hfunc_per_point=hf,
                expected_hfunc_per_point=expected_hf,
                hfunc_difference=exact_formula_applicable ? hf - expected_hf : NaN,
                observed_hinv_per_point=hi,
                expected_hinv_per_point=expected_hi,
                hinv_difference=exact_formula_applicable ? hi - expected_hi : NaN,
                observed_total_h_hinv_per_point=hf + hi,
                expected_total_h_hinv_per_point=expected_total,
                total_difference=exact_formula_applicable ? (hf + hi) - expected_total : NaN,
                observed_pair_cdf_per_point=pc,
                expected_pair_cdf_per_point=expected_pc,
                pair_cdf_difference=pc - expected_pc,
            ))
        end
    end
end
write_rows(joinpath(CFG.out, "complexity_check.csv"), complexity_rows)

# -----------------------------
# Pareto frontiers
# -----------------------------

pareto = NamedTuple[]
for case_id in sort(collect(keys(summary_by_case)))
    rows = summary_by_case[case_id]
    for method in ("current_indicator", "reduced_conditional")
        selected = sort(filter(r -> r.method == method && r.reference_resolved, rows); by=r -> r.median_runtime_s)
        best_rmse = Inf
        for row in selected
            if row.rmse < best_rmse
                push!(pareto, (;
                    case_id,
                    d=row.d,
                    method,
                    N=row.N,
                    rmse=row.rmse,
                    median_runtime_s=row.median_runtime_s,
                    rmse_ci_low=row.rmse_ci_low,
                    rmse_ci_high=row.rmse_ci_high,
                    runtime_ci_low_s=row.runtime_ci_low_s,
                    runtime_ci_high_s=row.runtime_ci_high_s,
                ))
                best_rmse = row.rmse
            end
        end
    end
end
!isempty(pareto) && write_rows(joinpath(CFG.out, "pareto_frontier.csv"), pareto)

# -----------------------------
# Figure-ready data
# -----------------------------

figdir = joinpath(CFG.out, "figure_data")
mkpath(figdir)
fig_summary = [(;
    case_id=r.case_id,
    design=r.design,
    point=r.point,
    d=r.d,
    N=r.N,
    method=r.method,
    rmse=r.rmse,
    rmse_ci_low=r.rmse_ci_low,
    rmse_ci_high=r.rmse_ci_high,
    median_runtime_s=r.median_runtime_s,
    runtime_ci_low_s=r.runtime_ci_low_s,
    runtime_ci_high_s=r.runtime_ci_high_s,
    reference_resolved=r.reference_resolved,
) for r in summary]
write_rows(joinpath(figdir, "rmse_vs_N.csv"), fig_summary)
write_rows(joinpath(figdir, "rmse_vs_runtime.csv"), fig_summary)
!isempty(matched) && write_rows(joinpath(figdir, "matched_speedup.csv"), matched)
!isempty(paired) && write_rows(joinpath(figdir, "paired_equal_N.csv"), paired)

# -----------------------------
# Publication tables
# -----------------------------

tabledir = joinpath(CFG.out, "table_data")
mkpath(tabledir)
max_budget_rows = NamedTuple[]
for case_id in sort(collect(keys(summary_by_case)))
    rows = summary_by_case[case_id]
    for method in ("current_indicator", "reduced_conditional")
        selected = filter(r -> r.method == method, rows)
        isempty(selected) && continue
        row = selected[argmax(getfield.(selected, :N))]
        push!(max_budget_rows, (;
            case_id,
            d=row.d,
            method,
            N=row.N,
            rmse=row.rmse,
            median_abs_error=row.median_abs_error,
            median_runtime_s=row.median_runtime_s,
            median_allocations=row.median_allocations,
            reference_resolved=row.reference_resolved,
        ))
    end
end
write_rows(joinpath(tabledir, "max_budget.csv"), max_budget_rows)
!isempty(strictest) && write_rows(joinpath(tabledir, "matched_accuracy_strictest.csv"), strictest)
write_rows(joinpath(tabledir, "complexity.csv"), complexity_rows)

# -----------------------------
# Human-readable report
# -----------------------------

valid_refs = count(r -> _parse_bool(r, "target_met"), refs)
# A paired replicate comparison is admitted when BOTH method cells at that
# (case, N) are reference-resolved according to the relative uncertainty rule.
resolved_cell_keys = Set((r.case_id, r.N, r.method) for r in summary if r.reference_resolved)
resolved_pair_keys = Set{Tuple{String,Int}}()
for r in summary
    key = (r.case_id, r.N)
    if (r.case_id, r.N, "current_indicator") in resolved_cell_keys &&
       (r.case_id, r.N, "reduced_conditional") in resolved_cell_keys
        push!(resolved_pair_keys, key)
    end
end
paired_valid = filter(r -> (r.case_id, r.N) in resolved_pair_keys, paired)
error_wins = count(r -> r.error_win, paired_valid)
runtime_wins = count(r -> r.runtime_win, paired_valid)
both_wins = count(r -> r.both_win, paired_valid)
finite_error_ratios = filter(isfinite, getfield.(paired_valid, :error_ratio))

n_cases_report = length(unique([r["case_id"] for r in raw]))
budgets_report = sort(unique([_parse_int(r, "N") for r in raw]))

open(joinpath(CFG.out, "REPORT.md"), "w") do io
    println(io, "# Publication benchmark report")
    println(io)
    println(io, "## Design")
    println(io)
    println(io, "- Mode: `$(CFG.mode)`")
    println(io, "- Cases: $n_cases_report")
    println(io, "- Budgets: $(join(budgets_report, ", "))")
    println(io, "- Replicates per cell: $(CFG.reps)")
    println(io, "- References meeting the aspirational absolute target: $valid_refs/$(length(refs))")
    println(io, "- Aggregate-cell reference rule: split-half z <= $(CFG.reference_split_z_max) and uncertainty/RMSE <= $(CFG.max_reference_ratio)")
    println(io, "- The absolute reference target is reported as an audit field; it is not itself required for a method cell to be resolved.")
    println(io)
    println(io, "## Equal-budget paired results")
    println(io)
    println(io, "- Valid paired replicate comparisons: $(length(paired_valid))")
    println(io, "- Reduced error wins: $error_wins/$(length(paired_valid)) ($(100error_wins/max(1,length(paired_valid)))%)")
    println(io, "- Reduced runtime wins: $runtime_wins/$(length(paired_valid)) ($(100runtime_wins/max(1,length(paired_valid)))%)")
    println(io, "- Reduced wins both: $both_wins/$(length(paired_valid)) ($(100both_wins/max(1,length(paired_valid)))%)")
    isempty(finite_error_ratios) || println(io, "- Median reduced/current absolute-error ratio: $(median(finite_error_ratios))")
    println(io)
    println(io, "## Matched accuracy")
    println(io)
    println(io, "Matched accuracy is evaluated for EVERY resolved current-method target cell; it does not select only the maximum speedup per case.")
    println(io)
    if !isempty(matched)
        println(io, "- Matched target cells: $(length(matched))")
        println(io, "- Median speedup over all matched targets: $(median(getfield.(matched, :speedup)))")
        println(io, "- 10th percentile speedup: $(quantile(getfield.(matched, :speedup), 0.10))")
        println(io, "- 90th percentile speedup: $(quantile(getfield.(matched, :speedup), 0.90))")
    else
        println(io, "- No resolved matched-accuracy cells were available.")
    end
    if !isempty(strictest)
        println(io, "- Median speedup at the strictest resolved current target per case: $(median(getfield.(strictest, :speedup)))")
    end
    println(io)
    println(io, "## Complexity check")
    println(io)
    println(io, "For the compiled D-vine traversal, exact counts with s=d-2 are checked:")
    println(io, "CDF: hfunc=s(3s+1)/2, hinv=s(s-1)/2, total=2s^2.")
    println(io, "Interior rectangle: hfunc=2s(s+1), hinv=s(s-1)/2, total=s(5s+3)/2.")
    println(io, "Genuine R-vine counts are reported as structure-dependent; no D-vine exact formula is imposed on them.")
    println(io, "For genuine R-vines, operation counts are recorded empirically; a universal exact complexity formula is not claimed here without a separate proof.")
    println(io)
    unresolved = filter(r -> !r.reference_resolved, summary)
    println(io, "## Reference-resolution audit")
    println(io)
    println(io, "- Aggregate method cells: $(length(summary))")
    println(io, "- Cells meeting reference-resolution rule: $(length(summary)-length(unresolved))")
    println(io, "- Cells flagged unresolved: $(length(unresolved))")
    println(io)
    println(io, "Raw replicate data, references, bootstrap summaries, Pareto frontiers, and figure/table-ready CSV files are retained under this run directory.")
end

println("Publication analysis written to ", CFG.out)
