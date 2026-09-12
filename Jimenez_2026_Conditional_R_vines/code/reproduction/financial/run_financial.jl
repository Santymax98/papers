using Copulas
using Dates
using Downloads
using Printf
using Random
using Statistics
using VineCopulas

include(joinpath(@__DIR__, "..", "..", "vendor", "VineCopulas", "benchmarks", "publication", "publication_common.jl"))

const OUTDIR = joinpath(@__DIR__, "rerun")
mkpath(OUTDIR)
const DATA_DIR = joinpath(@__DIR__, "..", "..", "data")
const URL = "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/10_Industry_Portfolios_daily_CSV.zip"
const ZIP_PATH = joinpath(DATA_DIR, "10_Industry_Portfolios_daily_CSV.zip")
const RAW_OUT = joinpath(OUTDIR, "05_kenneth_french_raw.csv")
const SUMMARY_OUT = joinpath(OUTDIR, "05_kenneth_french_summary.csv")
const REPORT_OUT = joinpath(OUTDIR, "05_kenneth_french_application.md")
const BASE_SEED = parse(Int, get(ENV, "KF_SEED", "20260905"))
const REPS = parse(Int, get(ENV, "KF_REPS", "30"))
const POWERS = parse.(Int, split(get(ENV, "KF_POWERS", "8,10,12"), ','))
const FIT_START = Date(get(ENV, "KF_START", "2015-01-01"))
const FIT_END = Date(get(ENV, "KF_END", "2024-12-31"))

function ensure_data!()
    mkpath(DATA_DIR)
    if !isfile(ZIP_PATH)
        Downloads.download(URL, ZIP_PATH)
    end
    return ZIP_PATH
end

function read_zip_text(path)
    return read(`unzip -p $path`, String)
end

function parse_industry_daily(text)
    lines = split(text, '\n')
    header_idx = findfirst(line -> startswith(strip(line), ",NoDur"), lines)
    header_idx === nothing && error("could not locate 10 Industry daily header")
    names = split(strip(lines[header_idx]), ',')[2:end]
    dates = Date[]
    data = Vector{Float64}[]
    for line in Iterators.drop(lines, header_idx)
        fields = split(strip(line), ',')
        length(fields) == 11 || continue
        occursin(r"^\d{8}$", fields[1]) || continue
        vals = try
            parse.(Float64, fields[2:end]) ./ 100
        catch
            continue
        end
        any(v -> v <= -0.99 || !isfinite(v), vals) && continue
        push!(dates, Date(fields[1], dateformat"yyyymmdd"))
        push!(data, vals)
    end
    X = Matrix{Float64}(undef, length(data), length(names))
    @inbounds for i in eachindex(data), j in eachindex(names)
        X[i, j] = data[i][j]
    end
    return names, dates, X
end

function pseudo_observations(X)
    n, d = size(X)
    U = Matrix{Float64}(undef, d, n)
    @inbounds for j in 1:d
        ord = sortperm(view(X, :, j); alg=MergeSort)
        ranks = zeros(Float64, n)
        i = 1
        while i <= n
            k = i
            while k < n && X[ord[k + 1], j] == X[ord[i], j]
                k += 1
            end
            avg = (i + k) / 2
            for t in i:k
                ranks[ord[t]] = avg
            end
            i = k + 1
        end
        for i in 1:n
            U[j, i] = ranks[i] / (n + 1)
        end
    end
    return U
end

function current_from_seed(vine, lower, upper, N, seed)
    Z = VineCopulas._reduced_owen_points(length(vine), N, UInt64(seed))
    return rectangle_indicator(inverse_rosenblatt(vine, Z), lower, upper)
end

function reduced_from_seed!(runner, N, seed)
    return reduced_compiled_estimate!(runner, N, seed)
end

function summarize(rows)
    groups = Dict{Tuple{String,Int,String},Vector{NamedTuple}}()
    for r in rows
        push!(get!(groups, (r.event, r.N, r.method), NamedTuple[]), r)
    end
    out = NamedTuple[]
    for key in sort(collect(keys(groups)))
        rr = groups[key]
        vals = getfield.(rr, :estimate)
        runt = getfield.(rr, :runtime_s)
        push!(out, (;
            event=key[1], N=key[2], method=key[3], n=length(rr),
            mean_estimate=mean(vals), sd_estimate=length(vals) > 1 ? std(vals) : NaN,
            rqmc_se=length(vals) > 1 ? std(vals) / sqrt(length(vals)) : NaN,
            median_runtime_s=median(runt),
            median_allocated_bytes=median(getfield.(rr, :bytes)),
            median_allocations=median(getfield.(rr, :allocations)),
            median_hfunc_calls=median(getfield.(rr, :hfunc_calls)),
            median_hinv_calls=median(getfield.(rr, :hinv_calls)),
            median_pair_cdf_calls=median(getfield.(rr, :pair_cdf_calls)),
        ))
    end
    return out
end

function main()
    zip = ensure_data!()
    names, dates, Xall = parse_industry_daily(read_zip_text(zip))
    keep = findall(d -> FIT_START <= d <= FIT_END, dates)
    length(keep) >= 1000 || error("insufficient Kenneth French observations in selected window")
    X = Xall[keep, :]
    kept_dates = dates[keep]
    U = pseudo_observations(X)

    fit_timed = @timed fit(
        RVineCopula, U;
        family_set=(GaussianCopula,),
        include_independence=true,
        allow_rotations=false,
        preselect=true,
        selection_criterion=:bic,
        tree_criterion=:tau,
        trace=false,
    )
    vine = fit_timed.value
    d = length(vine)
    lower_tail = fill(0.0, d)
    upper_tail = fill(0.10, d)
    central_lower = fill(0.25, d)
    central_upper = fill(0.75, d)
    events = (
        (event="joint_lower_10pct", lower=lower_tail, upper=upper_tail),
        (event="central_25_75_rectangle", lower=central_lower, upper=central_upper),
    )

    rows = NamedTuple[]
    for ev in events
        runner = prepare_reduced_runner((vine=vine, lower=ev.lower, upper=ev.upper))
        for pwr in POWERS
            N = 1 << pwr
            for rep in 1:REPS
                seed = BASE_SEED + 10_000 * pwr + rep
                current_t = @timed current_from_seed(vine, ev.lower, ev.upper, N, seed)
                push!(rows, (;
                    event=ev.event, d, N, rep, seed, method="current_indicator",
                    estimate=current_t.value, runtime_s=current_t.time,
                    bytes=current_t.bytes, allocations=Base.gc_alloc_count(current_t.gcstats),
                    hfunc_calls=NaN, hinv_calls=NaN, pair_cdf_calls=NaN,
                ))
                reduced_t = @timed reduced_from_seed!(runner, N, seed)
                push!(rows, (;
                    event=ev.event, d, N, rep, seed, method="reduced_conditional",
                    estimate=reduced_t.value.estimate, runtime_s=reduced_t.time,
                    bytes=reduced_t.bytes, allocations=Base.gc_alloc_count(reduced_t.gcstats),
                    hfunc_calls=reduced_t.value.hfunc_calls,
                    hinv_calls=reduced_t.value.hinv_calls,
                    pair_cdf_calls=reduced_t.value.pair_cdf_calls,
                ))
            end
        end
    end
    write_rows(RAW_OUT, rows)
    summary = summarize(rows)
    write_rows(SUMMARY_OUT, summary)

    open(REPORT_OUT, "w") do io
        println(io, "# Phase 5 Kenneth French financial application")
        println(io)
        println(io, "Status: COMPLETE illustrative application.")
        println(io)
        println(io, "## Data")
        println(io)
        println(io, "- Source: Kenneth R. French Data Library, `10_Industry_Portfolios_daily_CSV.zip`.")
        println(io, "- URL: `$URL`.")
        println(io, "- Cached file: `data/10_Industry_Portfolios_daily_CSV.zip`.")
        println(io, "- Sample window: $FIT_START to $FIT_END.")
        println(io, "- Observations used: $(length(keep)).")
        println(io, "- First/last retained date: $(first(kept_dates)) / $(last(kept_dates)).")
        println(io, "- Portfolios: $(join(names, ", ")).")
        println(io)
        println(io, "Returns were converted to continuous pseudo-observations by mid-ranks divided by n+1.")
        println(io)
        println(io, "## Model")
        println(io)
        println(io, "A full simplified Gaussian/independence R-vine was fitted by the package's sequential R-vine fitting routine. This is a copula-computation application, not an econometric claim that Gaussian pair-copulas are the best financial model.")
        println(io)
        println(io, "- Dimension: $d.")
        println(io, "- Fitting time: $(fit_timed.time) seconds.")
        println(io, "- Terminal pair type: $(typeof(VineCopulas._compile_reduced_probability(vine).terminal_copula)).")
        println(io)
        println(io, "## Events")
        println(io)
        println(io, "- `joint_lower_10pct`: `P(U_i <= 0.10, i=1,...,10)`.")
        println(io, "- `central_25_75_rectangle`: `P(0.25 < U_i <= 0.75, i=1,...,10)`.")
        println(io)
        println(io, "## Summary")
        println(io)
        println(io, "| event | N | method | mean estimate | RQMC SE | median runtime s |")
        println(io, "|---|---:|---|---:|---:|---:|")
        for r in summary
            println(io, "| $(r.event) | $(r.N) | $(r.method) | $(r.mean_estimate) | $(r.rqmc_se) | $(r.median_runtime_s) |")
        end
        println(io)
        println(io, "The same seeds and Sobol/Owen construction are used for paired current/reduced evaluations. No public `cdf` API was changed.")
    end
    println("Wrote Kenneth French application outputs to $OUTDIR")
end

main()
