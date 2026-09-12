include("publication_common.jl")

const CFG = publication_config()
path = joinpath(CFG.out, "external_rvinecopulib.csv")
rows = read_table(path)
isempty(rows) && error("missing/empty external_rvinecopulib.csv; run rvinecopulib_external.R first")

groups = Dict{Tuple{String,Int},Vector{Dict{String,String}}}()
for row in rows
    key = (row["case_id"], _parse_int(row, "N"))
    push!(get!(groups, key, Dict{String,String}[]), row)
end

out = NamedTuple[]
for key in sort!(collect(keys(groups)); by=string)
    rs = groups[key]
    errors = [_parse_float(r, "abs_error") for r in rs]
    runtimes = [_parse_float(r, "runtime_s") for r in rs]
    estimates = [_parse_float(r, "estimate") for r in rs]
    push!(out, (;
        case_id=key[1],
        d=_parse_int(first(rs), "d"),
        N=key[2],
        n=length(rs),
        mean_estimate=mean(estimates),
        estimate_sd=length(estimates) > 1 ? std(estimates) : 0.0,
        median_abs_error=median(errors),
        rmse=sqrt(mean(abs2, errors)),
        median_runtime_s=median(runtimes),
    ))
end

write_rows(joinpath(CFG.out, "external_rvinecopulib_summary.csv"), out)
println("External rvinecopulib analysis written to ", CFG.out)
