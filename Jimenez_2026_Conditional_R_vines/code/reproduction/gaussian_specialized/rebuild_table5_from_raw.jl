using Printf
using Statistics

const HERE = @__DIR__
const PACKAGE_ROOT = normpath(joinpath(HERE, "..", "..", ".."))
const PUBLICATION_COMMON = joinpath(PACKAGE_ROOT, "code", "benchmarks", "publication", "publication_common.jl")
include(PUBLICATION_COMMON)

const METHOD_RANK = Dict(
    "current_indicator" => 1,
    "reduced_conditional" => 2,
    "additionaldistributions" => 3,
    "mvtnorm_genzbretz" => 4,
)
const METHOD_LABEL = Dict(
    "current_indicator" => "Indicator",
    "reduced_conditional" => "Reduced",
    "additionaldistributions" => "AdditionalDist.",
    "mvtnorm_genzbretz" => "mvtnorm Genz--Bretz",
)

raw = vcat(
    read_table(joinpath(HERE, "raw_results.csv")),
    read_table(joinpath(HERE, "mvtnorm_raw_results.csv")),
)
isempty(raw) && error("Gaussian raw results are missing")

groups = Dict{Tuple{String,String},Vector{Dict{String,String}}}()
for row in raw
    push!(get!(groups, (row["case_id"], row["method"]), Dict{String,String}[]), row)
end

summary = NamedTuple[]
for key in sort(collect(keys(groups)); by = key -> (
    _parse_int(first(groups[key]), "d"), key[1], METHOD_RANK[key[2]],
))
    rows = groups[key]
    reference = _parse_float(first(rows), "reference")
    estimates = [_parse_float(row, "estimate") for row in rows]
    runtimes = [_parse_float(row, "runtime_s") for row in rows]
    push!(summary, (;
        case_id=key[1],
        d=_parse_int(first(rows), "d"),
        event=first(rows)["event"],
        method=key[2],
        rmse=sqrt(mean((estimates .- reference) .^ 2)),
        median_runtime_s=median(runtimes),
    ))
end

path = get(ENV, "OUT_TABLE", joinpath(PACKAGE_ROOT, "output", "gaussian_specialized_comparison.tex"))
mkpath(dirname(path))
open(path, "w") do io
    println(io, "\\begin{table}[H]")
    println(io, "\\centering\\small")
    println(io, "\\begin{tabular}{lrllrr}")
    println(io, "\\toprule")
    println(io, "Case & \$d\$ & Event & Method & RMSE & Median time (ms) \\\\")
    println(io, "\\midrule")
    previous_case = ""
    for row in summary
        !isempty(previous_case) && previous_case != row.case_id && println(io, "\\addlinespace")
        case_label = replace(row.case_id, "gaussian_" => "", "_d$(row.d)" => "", "_" => " ")
        @printf(io, "%s & %d & %s & %s & %.3e & %.2f",
            case_label, row.d, row.event, METHOD_LABEL[row.method], row.rmse,
            1000 * row.median_runtime_s)
        println(io, " \\\\")
        previous_case = row.case_id
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
    println(io, "\\caption{Gaussian-specific comparison at nominal budget/max-points 4096 and 30 independent randomizations. RMSE is measured against the frozen deterministic cubature reference. Reduced and Indicator use scrambled Sobol nets; \\texttt{AdditionalDistributions.jl} 0.3.0 uses 10 random shifts with MVSORT conditioning and a CBC rank-1 lattice; \\texttt{mvtnorm} 1.3-3 uses \\texttt{GenzBretz}. Both specialized calls request absolute and relative tolerances of \$10^{-6}\$. JIT and one-time construction costs are excluded.}")
    println(io, "\\label{tab:s-gaussian-specialized}")
    println(io, "\\end{table}")
end
println("Table 5 reconstructed from ", length(raw), " raw rows: ", path)

