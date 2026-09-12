include("publication_common.jl")

const CFG = publication_config()
const CASES = publication_cases(CFG)
const SAMPLES = _env_int("REGULARITY_SAMPLES", CFG.mode == "final" ? 512 : 64)
const STEP = _env_float("REGULARITY_STEP", 2e-6)
const REG_SEED = _env_int("REGULARITY_SEED", CFG.base_seed + 400_000)

mkpath(CFG.out)
raw_path = joinpath(CFG.out, "regularity_raw.csv")
isfile(raw_path) && rm(raw_path; force=true)

function transported_component!(runner::PublicationReducedRunner, z, variable)
    VineCopulas._truncated_inverse_rosenblatt!(
        runner.conditioning_values,
        runner.workspace,
        runner.plan,
        runner.lower,
        runner.upper,
        z,
        runner.order,
        runner.slots,
    )
    return runner.conditioning_values[variable]
end

summary_rows = NamedTuple[]
println("Regularity diagnostic: samples=$SAMPLES step=$STEP")

for (case_idx, case) in enumerate(CASES)
    runner = prepare_reduced_runner(case)
    s = length(case.vine) - 2
    Z = VineCopulas._reduced_owen_points(s, SAMPLES, UInt64(REG_SEED + 10_000case_idx))
    case_kappa = Float64[]
    stage_kappa = [Float64[] for _ in 1:s]
    zp = Vector{Float64}(undef, s)
    zm = Vector{Float64}(undef, s)

    for col in axes(Z, 2)
        z = @view Z[:, col]
        for stage in 1:s
            h = min(STEP, 0.45 * min(z[stage], 1 - z[stage]))
            h > 0 || continue
            copyto!(zp, z)
            copyto!(zm, z)
            zp[stage] += h
            zm[stage] -= h
            variable = runner.order[stage]
            vp = transported_component!(runner, zp, variable)
            vm = transported_component!(runner, zm, variable)
            kappa = abs(vp - vm) / (2h)
            isfinite(kappa) || continue
            push!(case_kappa, kappa)
            push!(stage_kappa[stage], kappa)
            append_row(raw_path, (;
                case_id=case.id,
                d=length(case.vine),
                sample=col,
                stage,
                variable,
                z_stage=z[stage],
                kappa,
                finite_difference_step=h,
            ))
        end
    end

    for stage in 1:s
        vals = stage_kappa[stage]
        isempty(vals) && continue
        push!(summary_rows, (;
            case_id=case.id,
            d=length(case.vine),
            stage,
            variable=runner.order[stage],
            n=length(vals),
            median_kappa=median(vals),
            p90_kappa=quantile(vals, 0.90),
            p99_kappa=quantile(vals, 0.99),
            max_kappa=maximum(vals),
        ))
    end
    @printf("[%d/%d] %-34s max kappa=%.5g\n", case_idx, length(CASES), case.id, maximum(case_kappa))
end

write_rows(joinpath(CFG.out, "regularity_by_stage.csv"), summary_rows)
println("Regularity diagnostics written to ", CFG.out)
