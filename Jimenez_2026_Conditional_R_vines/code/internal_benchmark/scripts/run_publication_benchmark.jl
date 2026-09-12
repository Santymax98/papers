include("publication_common.jl")

const CFG = publication_config()
const CASES = publication_cases(CFG)

mkpath(CFG.out)
capture_environment!(CFG, CASES)
refs = aggregate_reference_files!(CFG, CASES)
ref_by_case = Dict(r.case_id => r for r in refs)
rawdir = joinpath(CFG.out, "raw")
mkpath(rawdir)

function completed_keys(path)
    keys = Set{Tuple{Int,Int,String}}()
    for row in read_table(path)
        push!(keys, (_parse_int(row, "N"), _parse_int(row, "rep"), row["method"]))
    end
    return keys
end

function collect_raw_files!(cfg, cases)
    rows = Dict{String,String}[]
    for case in cases
        append!(rows, read_table(joinpath(cfg.out, "raw", case.id * ".csv")))
    end
    isempty(rows) && error("no raw publication rows found")
    # Recreate typed rows in a stable column order rather than concatenating
    # arbitrary text files.
    typed = NamedTuple[]
    for r in rows
        push!(typed, (;
            case_id=r["case_id"],
            design=r["design"],
            model=r["model"],
            point=r["point"],
            d=_parse_int(r, "d"),
            conditioning_dimension=_parse_int(r, "conditioning_dimension"),
            power=_parse_int(r, "power"),
            N=_parse_int(r, "N"),
            rep=_parse_int(r, "rep"),
            seed=_parse_int(r, "seed"),
            method=r["method"],
            execution_rank=_parse_int(r, "execution_rank"),
            estimate=_parse_float(r, "estimate"),
            reference=_parse_float(r, "reference"),
            reference_uncertainty=_parse_float(r, "reference_uncertainty"),
            reference_target_met=_parse_bool(r, "reference_target_met"),
            abs_error=_parse_float(r, "abs_error"),
            runtime_s=_parse_float(r, "runtime_s"),
            bytes=_parse_float(r, "bytes"),
            allocations=_parse_float(r, "allocations"),
            hfunc_calls=_parse_float(r, "hfunc_calls"),
            hinv_calls=_parse_float(r, "hinv_calls"),
            pair_cdf_calls=_parse_float(r, "pair_cdf_calls"),
        ))
    end
    sort!(typed; by=r -> (r.d, r.case_id, r.N, r.rep, r.method))
    write_rows(joinpath(cfg.out, "raw_results.csv"), typed)
    return typed
end

println("Publication benchmark")
println("  mode   = ", CFG.mode)
println("  cases  = ", length(CASES))
println("  powers = ", join(CFG.powers, ','))
println("  reps   = ", CFG.reps)
println("  output = ", CFG.out)
println("  Julia threads = ", Threads.nthreads())
Threads.nthreads() == 1 || @warn "Publication timing is intended for julia -t 1"

case_costs = NamedTuple[]

for (case_idx, case) in enumerate(CASES)
    ref = ref_by_case[case.id]
    rawpath = joinpath(rawdir, case.id * ".csv")
    done = CFG.resume ? completed_keys(rawpath) : Set{Tuple{Int,Int,String}}()
    !CFG.resume && isfile(rawpath) && rm(rawpath; force=true)

    # Warm JIT separately, then measure algorithmic plan construction after JIT.
    prepare_reduced_runner(case)
    plan_timed = @timed prepare_reduced_runner(case)
    runner = plan_timed.value
    push!(case_costs, (;
        case_id=case.id,
        d=length(case.vine),
        plan_build_runtime_s=plan_timed.time,
        plan_build_bytes=plan_timed.bytes,
        plan_build_allocations=_gc_allocs(plan_timed),
        order=join(runner.order, '-'),
    ))

    # Warm both estimators. These calls are deliberately not recorded.
    current_estimate(case.vine, case.lower, case.upper, 64, CFG.base_seed)
    reduced_compiled_estimate!(runner, 64, CFG.base_seed)

    @printf("[%d/%d] %s\n", case_idx, length(CASES), case.id)
    for power in CFG.powers
        N = 1 << power
        GC.gc()
        for rep in 1:CFG.reps
            seed = CFG.base_seed + rep - 1
            # Alternate method order to reduce systematic thermal/cache bias.
            method_order = isodd(rep) ? ("current_indicator", "reduced_conditional") :
                                        ("reduced_conditional", "current_indicator")
            for (rank, method) in enumerate(method_order)
                key = (N, rep, method)
                key in done && continue

                if method == "current_indicator"
                    timed = @timed current_estimate(case.vine, case.lower, case.upper, N, seed)
                    estimate = timed.value
                    hfunc_calls = NaN
                    hinv_calls = NaN
                    pair_cdf_calls = NaN
                else
                    timed = @timed reduced_compiled_estimate!(runner, N, seed)
                    result = timed.value
                    estimate = result.estimate
                    hfunc_calls = result.hfunc_calls
                    hinv_calls = result.hinv_calls
                    pair_cdf_calls = result.pair_cdf_calls
                end

                row = (;
                    case_id=case.id,
                    design=case.design,
                    model=case.model,
                    point=case.point,
                    d=length(case.vine),
                    conditioning_dimension=length(case.vine)-2,
                    power,
                    N,
                    rep,
                    seed,
                    method,
                    execution_rank=rank,
                    estimate,
                    reference=ref.estimate,
                    reference_uncertainty=ref.uncertainty,
                    reference_target_met=ref.target_met,
                    abs_error=abs(estimate - ref.estimate),
                    runtime_s=timed.time,
                    bytes=timed.bytes,
                    allocations=_gc_allocs(timed),
                    hfunc_calls,
                    hinv_calls,
                    pair_cdf_calls,
                )
                append_row(rawpath, row)
                push!(done, key)
            end
        end
        @printf("    N=%d complete\n", N)
        flush(stdout)
    end
end

write_rows(joinpath(CFG.out, "plan_costs.csv"), case_costs)
collect_raw_files!(CFG, CASES)
println("Raw publication benchmark complete: ", CFG.out)
