include("publication_common.jl")

const CFG = publication_config()
const CASES = publication_cases(CFG)

mkpath(CFG.out)
mkpath(reference_dir(CFG))
capture_environment!(CFG, CASES)

function save_reference!(case, data::Dict{String,Any})
    path = reference_toml(CFG, case.id)
    open(path, "w") do io
        TOML.print(io, data)
    end
    return path
end

function compute_hcubature_reference(case)
    f, s = reference_integrand(case)
    timed = @timed Copulas.HCubature.hcubature(
        f,
        zeros(s),
        ones(s);
        rtol=CFG.reference_hcubature_rtol,
        atol=CFG.reference_hcubature_atol,
        maxevals=CFG.reference_maxevals,
    )
    estimate, error_estimate = timed.value
    target_met = error_estimate <= CFG.reference_target_se
    return Dict{String,Any}(
        "case_id" => case.id,
        "d" => length(case.vine),
        "conditioning_dimension" => s,
        "method" => "deterministic_hcubature_reduced_integrand",
        "estimate" => Float64(estimate),
        "uncertainty" => Float64(error_estimate),
        "uncertainty_kind" => "hcubature_error_estimate",
        "N" => 0,
        "R" => 0,
        "seed" => 0,
        "target_met" => target_met,
        "split_half_z" => 0.0,
        "runtime_s" => timed.time,
        "bytes" => timed.bytes,
        "allocations" => _gc_allocs(timed),
    )
end

function compute_rqmc_reference(case, case_idx::Int)
    history_path = joinpath(reference_dir(CFG), case.id * "_history.csv")
    replicates_path = joinpath(reference_dir(CFG), case.id * "_replicates.csv")
    isfile(history_path) && rm(history_path; force=true)
    isfile(replicates_path) && rm(replicates_path; force=true)

    final = nothing
    for power in CFG.reference_powers
        N = 1 << power
        seed = CFG.reference_seed + 1_000_000 * case_idx + 10_000 * power
        timed = @timed VineCopulas._rectprob_reduced_rqmc(
            case.vine,
            case.lower,
            case.upper;
            N=N,
            R=CFG.reference_r,
            seed=seed,
            order_strategy=:default,
            randomized=true,
        )
        result = timed.value
        reps = result.replicates
        m1, m2, se1, se2, split_z = split_half_diagnostic(reps)
        target_met = isfinite(result.stderr) &&
                     result.stderr <= CFG.reference_target_se &&
                     split_z <= CFG.reference_split_z_max

        append_row(history_path, (;
            case_id=case.id,
            d=length(case.vine),
            power,
            N,
            R=CFG.reference_r,
            seed,
            estimate=result.estimate,
            stderr=result.stderr,
            split_mean_1=m1,
            split_mean_2=m2,
            split_se_1=se1,
            split_se_2=se2,
            split_z,
            target_met,
            runtime_s=timed.time,
            bytes=timed.bytes,
            allocations=_gc_allocs(timed),
        ))
        for (rep, estimate) in enumerate(reps)
            append_row(replicates_path, (;
                case_id=case.id,
                power,
                N,
                rep,
                seed=seed + rep - 1,
                estimate,
            ))
        end

        final = Dict{String,Any}(
            "case_id" => case.id,
            "d" => length(case.vine),
            "conditioning_dimension" => length(case.vine) - 2,
            "method" => "independent_reduced_rqmc_reference",
            "estimate" => result.estimate,
            "uncertainty" => result.stderr,
            "uncertainty_kind" => "between_scramble_standard_error",
            "N" => N,
            "R" => CFG.reference_r,
            "seed" => seed,
            "target_met" => target_met,
            "split_half_z" => split_z,
            "split_mean_1" => m1,
            "split_mean_2" => m2,
            "runtime_s" => timed.time,
            "bytes" => timed.bytes,
            "allocations" => _gc_allocs(timed),
        )
        target_met && break
    end
    return final
end

println("Publication reference generation")
println("  mode       = ", CFG.mode)
println("  cases      = ", length(CASES))
println("  target SE  = ", CFG.reference_target_se)
println("  output     = ", CFG.out)

for (idx, case) in enumerate(CASES)
    path = reference_toml(CFG, case.id)
    if isfile(path) && CFG.resume && !CFG.force_references
        existing = TOML.parsefile(path)
        if get(existing, "target_met", false)
            @printf("[%d/%d] %-34s reuse target-met reference\n", idx, length(CASES), case.id)
            continue
        end
    end

    p = length(case.vine)
    # The reduced integral has dimension d-2; deterministic cubature is kept
    # only through dimension three to retain a genuinely independent numerical
    # reference in low dimensions.
    data = if p <= 5
        compute_hcubature_reference(case)
    else
        compute_rqmc_reference(case, idx)
    end
    save_reference!(case, data)
    @printf(
        "[%d/%d] %-34s estimate=%.12g uncertainty=%.3e target=%s\n",
        idx,
        length(CASES),
        case.id,
        data["estimate"],
        data["uncertainty"],
        string(data["target_met"]),
    )
    flush(stdout)
end

aggregate_reference_files!(CFG, CASES)
println("References written under ", reference_dir(CFG))
