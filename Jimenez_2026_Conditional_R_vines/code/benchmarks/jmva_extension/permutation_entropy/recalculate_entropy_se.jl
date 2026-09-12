using Statistics
using VineCopulas

include(joinpath(@__DIR__, "..", "..", "publication", "publication_common.jl"))

const RESULTS = get(
    ENV,
    "RESULTS",
    joinpath(@__DIR__, "..", "..", "..", "reports", "jmva_extension",
             "permutation_entropy", "final"),
)
const RAW_FILE = joinpath(RESULTS, "ordinal_probabilities.csv")
const SUMMARY_FILE = joinpath(RESULTS, "permutation_entropy_summary.csv")
const AUDIT_FILE = joinpath(RESULTS, "permutation_entropy_se_correction.csv")

raw_rows = read_table(RAW_FILE)
summary_rows = read_table(SUMMARY_FILE)
isempty(raw_rows) && error("no ordinal-probability rows in $RAW_FILE")
isempty(summary_rows) && error("no entropy-summary rows in $SUMMARY_FILE")

groups = Dict{Tuple{String,String},Vector{Dict{String,String}}}()
for row in raw_rows
    push!(get!(groups, (row["industry"], row["transform"]), Dict{String,String}[]), row)
end

audit_rows = NamedTuple[]
for row in summary_rows
    key = (row["industry"], row["transform"])
    pattern_rows = sort(groups[key]; by=r -> parse(Int, r["pattern_index"]))
    probabilities = [parse(Float64, r["model_probability"]) for r in pattern_rows]
    standard_errors = [parse(Float64, r["numerical_se"]) for r in pattern_rows]
    entropy = VineCopulas._normalized_permutation_entropy_with_gradient(probabilities)
    corrected_se = sqrt(sum(abs2, entropy.gradient .* standard_errors))
    old_se = parse(Float64, row["permutation_entropy_numerical_se"])

    isapprox(
        entropy.entropy,
        parse(Float64, row["permutation_entropy_model"]);
        atol=2e-15,
        rtol=2e-15,
    ) || error("entropy point estimate changed for $(key)")
    row["permutation_entropy_numerical_se"] = repr(corrected_se)
    push!(audit_rows, (;
        industry=key[1],
        transform=key[2],
        old_diagonal_se=old_se,
        corrected_diagonal_se=corrected_se,
        corrected_over_old=corrected_se / old_se,
        radial_derivative=sum(entropy.gradient .* probabilities),
    ))
end

header = _csv_parse_line(first(readlines(SUMMARY_FILE)))
open(SUMMARY_FILE, "w") do io
    println(io, join((_csv_escape(name) for name in header), ','))
    for row in summary_rows
        println(io, join((_csv_escape(row[name]) for name in header), ','))
    end
end
write_rows(AUDIT_FILE, audit_rows)

println("Corrected diagonal delta-method entropy SEs: ", SUMMARY_FILE)
println("Audit table: ", AUDIT_FILE)
