include("publication_common.jl")

const CFG = publication_config()
const CASES = publication_cases(CFG)
const ORDER_N = _env_int("ORDER_N", 1 << 10)
const ORDER_REPS = _env_int("ORDER_REPS", CFG.mode == "final" ? 30 : 4)
const ORDER_SEED = _env_int("ORDER_SEED", CFG.base_seed + 600_000)

rows = NamedTuple[]
changed = 0

for case in CASES
    widths = case.upper .- case.lower
    default_runner = prepare_reduced_runner(case; order_strategy=:default)
    small_runner = prepare_reduced_runner(case; order_strategy=:smallest_width_first)
    default_runner.order == small_runner.order && continue
    global changed += 1
    ref = load_reference(CFG, case.id)
    reference = Float64(ref["estimate"])

    for rep in 1:ORDER_REPS
        seed = ORDER_SEED + rep - 1
        for (strategy, runner) in (("default", default_runner), ("smallest_width_first", small_runner))
            timed = @timed reduced_compiled_estimate!(runner, ORDER_N, seed)
            result = timed.value
            push!(rows, (;
                case_id=case.id,
                d=length(case.vine),
                strategy,
                order=join(runner.order, '-'),
                q1=widths[first(runner.order)],
                N=ORDER_N,
                rep,
                seed,
                estimate=result.estimate,
                abs_error=abs(result.estimate - reference),
                runtime_s=timed.time,
                bytes=timed.bytes,
                allocations=_gc_allocs(timed),
                hfunc_calls=result.hfunc_calls,
                hinv_calls=result.hinv_calls,
            ))
        end
    end
end

if isempty(rows)
    println("No publication case changed order under :smallest_width_first; no CSV written.")
else
    write_rows(joinpath(CFG.out, "sampling_order_study.csv"), rows)
    println("Sampling-order study written for $changed cases to ", CFG.out)
end
