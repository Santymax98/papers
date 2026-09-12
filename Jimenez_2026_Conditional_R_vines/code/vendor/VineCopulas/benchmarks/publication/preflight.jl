include("publication_common.jl")

const CFG = publication_config()
const CASES = publication_cases(CFG)
const PREFLIGHT_POINTS = _env_int("PREFLIGHT_POINTS", 16)

println("Publication preflight")
println("  git commit = ", _shell_output(["git", "rev-parse", "HEAD"]))
println("  git branch = ", _shell_output(["git", "branch", "--show-current"]))
println("  git status = ", repr(_shell_output(["git", "status", "--porcelain"])))
println("  Julia      = ", VERSION)
println("  threads    = ", Threads.nthreads())
Threads.nthreads() == 1 || @warn "Use julia -t 1 for publication timing"

for (idx, case) in enumerate(CASES)
    runner = prepare_reduced_runner(case)
    s = length(case.vine) - 2
    Z = VineCopulas._reduced_owen_points(s, PREFLIGHT_POINTS, UInt64(CFG.base_seed + 800_000 + idx))
    max_violation = 0.0
    for col in axes(Z, 2)
        G, q1 = VineCopulas._reduced_probability_integrand!(
            runner.workspace,
            runner.plan,
            runner.lower,
            runner.upper,
            view(Z, :, col),
            runner.order,
            runner.slots,
            runner.conditioning_values,
        )
        max_violation = max(max_violation, max(-G, G - q1, 0.0))
    end
    max_violation <= 4096eps(Float64) || error("integrand bound violation in $(case.id)")

    # One-point call-count audit.  Do not request an Owen-scrambled net with
    # N=1: QuasiMonteCarlo's nested Owen implementation has zero scrambling
    # bits in that degenerate case.  Call the compiled integrand directly at
    # one deterministic interior point instead.
    reset_counts!(runner)
    z_audit = fill(0.5, s)
    G_audit, q1_audit = VineCopulas._reduced_probability_integrand!(
        runner.workspace,
        runner.plan,
        runner.lower,
        runner.upper,
        z_audit,
        runner.order,
        runner.slots,
        runner.conditioning_values,
    )
    (0.0 <= G_audit <= q1_audit <= 1.0) || error(
        "one-point integrand audit failed for $(case.id)",
    )
    result = (;
        hfunc_calls=runner.workspace.counts.hfunc,
        hinv_calls=runner.workspace.counts.hinv,
        pair_cdf_calls=runner.workspace.counts.pair_cdf,
    )
    expected_pc = case.point == "rectangle" ? 4 : 1
    result.pair_cdf_calls == expected_pc || error(
        "pair-CDF count mismatch for $(case.id): $(result.pair_cdf_calls) != $expected_pc",
    )
    if case.model == "D"
        # Exact counts for this compiled D-vine traversal. General R-vines can
        # reuse states differently and may require fewer calls.
        expected_hi = s * (s - 1) ÷ 2
        expected_hf = case.point == "rectangle" ? 2s * (s + 1) : s * (3s + 1) ÷ 2
        result.hinv_calls == expected_hi || error(
            "D-vine hinv count mismatch for $(case.id): $(result.hinv_calls) != $expected_hi",
        )
        result.hfunc_calls == expected_hf || error(
            "D-vine hfunc count mismatch for $(case.id): $(result.hfunc_calls) != $expected_hf",
        )
    end
    @printf("[%d/%d] %-34s bound OK; calls=(%d,%d,%d)\n",
            idx, length(CASES), case.id,
            result.hfunc_calls, result.hinv_calls, result.pair_cdf_calls)
end


scramble_rows = NamedTuple[]
for d in sort(unique(length(c.vine) for c in CASES))
    sdim = d - 2
    for rep in 1:4
        seed = CFG.base_seed + 950_000 + rep - 1
        full = VineCopulas._reduced_owen_points(d, 64, UInt64(seed))
        reduced = VineCopulas._reduced_owen_points(sdim, 64, UInt64(seed))
        max_difference = maximum(abs.(@view(full[1:sdim, :]) .- reduced))
        push!(scramble_rows, (; d, conditioning_dimension=sdim, rep, seed, max_difference,
                                bitwise_equal=max_difference == 0.0))
    end
end
write_rows(joinpath(CFG.out, "scramble_audit.csv"), scramble_rows)
all(r -> r.bitwise_equal, scramble_rows) || error("paired scramble audit failed")
println("Paired scramble audit passed; first d-2 axes are bitwise identical.")

println("Preflight passed for all publication cases.")
