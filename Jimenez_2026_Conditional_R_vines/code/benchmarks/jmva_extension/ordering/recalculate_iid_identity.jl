include(joinpath(@__DIR__, "..", "..", "publication", "publication_common.jl"))

const RESULTS = get(
    ENV,
    "RESULTS",
    joinpath(@__DIR__, "..", "..", "..", "reports", "jmva_extension",
             "ordering", "final"),
)
const SUMMARY_FILE = joinpath(RESULTS, "ordering_summary.csv")
const IID_N = parse(Int, get(ENV, "IID_N", "1024"))

rows = read_table(SUMMARY_FILE)
isempty(rows) && error("no ordering summary rows in $SUMMARY_FILE")
for row in rows
    V = parse(Float64, row["V_pi"])
    rmse = parse(Float64, row["iid_rmse"])
    row["iid_N_rmse2_over_V"] = repr(IID_N * rmse^2 / V)
end

old_header = _csv_parse_line(first(readlines(SUMMARY_FILE)))
header = if "iid_N_rmse2_over_V" in old_header
    old_header
else
    index = findfirst(==("iid_rmse"), old_header)
    vcat(old_header[1:index], "iid_N_rmse2_over_V", old_header[(index + 1):end])
end
open(SUMMARY_FILE, "w") do io
    println(io, join((_csv_escape(name) for name in header), ','))
    for row in rows
        println(io, join((_csv_escape(row[name]) for name in header), ','))
    end
end

println("Added N*RMSE^2/V diagnostic: ", SUMMARY_FILE)
