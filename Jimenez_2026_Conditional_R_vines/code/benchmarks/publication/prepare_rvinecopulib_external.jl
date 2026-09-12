include("publication_common.jl")

const CFG = publication_config()
const CASES = publication_cases(CFG)
const EXT_DIR = joinpath(CFG.out, "external_rvinecopulib")
mkpath(EXT_DIR)

refs = read_table(joinpath(CFG.out, "references.csv"))
isempty(refs) && error("references.csv is missing/empty; generate publication references first")
ref_by_case = Dict(r["case_id"] => r for r in refs)

# External pvinecop comparison is intentionally CDF-only.  General rectangles
# would require 2^d inclusion-exclusion calls and would not be a fair baseline.
cases = filter(c -> all(iszero, c.lower), CASES)
isempty(cases) && error("no CDF publication cases found")

# Build deterministic audit points from *interior Rosenblatt coordinates*.
# Arbitrary coordinate-wise points can drive high-tree conditional states
# extremely close to 0/1 in a long heterogeneous vine even when every raw
# coordinate is visually interior.  That can expose unrelated numerical
# boundary issues in density evaluation.  Starting from central latent
# Rosenblatt coordinates keeps the conditional states used to generate the
# audit observations away from those artificial extremes while still yielding
# nontrivial, model-specific points in the original u-space.
function audit_latent(d::Int, which::Int)
    if which == 1
        return fill(0.5, d)
    elseif which == 2
        return [0.5 + 0.08 * sin(2pi * (j - 1) / d) for j in 1:d]
    end
    error("unknown audit point")
end

function audit_point(vine, d::Int, which::Int)
    z = audit_latent(d, which)
    x = inverse_rosenblatt(vine, z)
    all(isfinite, x) || error("non-finite inverse-Rosenblatt audit point")
    all(t -> 0.0 < t < 1.0, x) || error("boundary inverse-Rosenblatt audit point")
    return x, z
end

case_rows = NamedTuple[]
audit_rows = NamedTuple[]
for case in cases
    ref = get(ref_by_case, case.id, nothing)
    ref === nothing && error("missing reference for $(case.id)")
    d = length(case.vine)
    push!(case_rows, (;
        case_id=case.id,
        design=case.design,
        model=case.model,
        point=case.point,
        d,
        upper=_vector_field(case.upper),
        reference=_parse_float(ref, "estimate"),
        reference_uncertainty=_parse_float(ref, "uncertainty"),
        reference_target_met=_parse_bool(ref, "target_met"),
        reference_split_half_z=_parse_float(ref, "split_half_z"),
    ))

    # Deterministic model-parity audit.  CDF estimators are stochastic/QMC,
    # therefore density/log-density is the clean way to establish that Julia
    # and rvinecopulib encode the same structure, families, and parameters.
    for which in 1:2
        x, z = audit_point(case.vine, d, which)
        lp = logpdf(case.vine, x)
        isfinite(lp) || error("non-finite Julia logpdf audit value for $(case.id) at Rosenblatt audit point $which")
        push!(audit_rows, (;
            case_id=case.id,
            d,
            audit_point=which,
            latent_z=_vector_field(z),
            u=_vector_field(x),
            julia_logpdf=lp,
        ))
    end
end

write_rows(joinpath(EXT_DIR, "cases.csv"), case_rows)
write_rows(joinpath(EXT_DIR, "julia_density_audit.csv"), audit_rows)

println("Prepared rvinecopulib external benchmark")
println("  CDF cases = ", length(cases))
println("  audit rows = ", length(audit_rows))
println("  output     = ", EXT_DIR)
