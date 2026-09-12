using Printf
using Random
using Statistics
using TOML
using VineCopulas

include(joinpath(@__DIR__, "..", "..", "publication", "publication_common.jl"))

const OUT = get(ENV, "OUT", joinpath(@__DIR__, "results"))
const IID_VARIANCE_N = parse(Int, get(ENV, "IID_VARIANCE_N", string(1 << 18)))
const IID_N = parse(Int, get(ENV, "IID_N", string(1 << 10)))
const IID_REPS = parse(Int, get(ENV, "IID_REPS", "100"))
const RQMC_POWERS = parse.(Int, split(get(ENV, "RQMC_POWERS", "8,10,12,14"), ','))
const RQMC_REPS = parse(Int, get(ENV, "RQMC_REPS", "30"))
const ROUGHNESS_N = parse(Int, get(ENV, "ROUGHNESS_N", "256"))
const BASE_SEED = parse(Int, get(ENV, "BASE_SEED", "730001"))

mkpath(OUT)

function cases()
    return [
        (;
            id="gaussian_rectangle_d5",
            vine=publication_homogeneous_dvine(5, 0.45),
            lower=fill(0.08, 5),
            upper=collect(range(0.48, 0.88; length=5)),
            interior=true,
        ),
        (;
            id="gaussian_strong_d5",
            vine=publication_homogeneous_dvine(5, 0.82),
            lower=zeros(5),
            upper=fill(0.72, 5),
            interior=false,
        ),
        (;
            id="genuine_rvine_d5",
            vine=publication_genuine_rvine5(),
            lower=zeros(5),
            upper=[0.36, 0.61, 0.47, 0.73, 0.58],
            interior=false,
        ),
    ]
end

function reference(case_id)
    file = joinpath(@__DIR__, "frozen_references", "$case_id.toml")
    data = TOML.parsefile(file)
    return Float64(data["estimate"]), Float64(data["uncertainty"])
end

order_string(order) = join(order, '-')

function integrand_samples!(out, case, plan, order, slots, Z)
    ws = VineCopulas._reduced_workspace(plan)
    values = fill(NaN, length(case.vine))
    @inbounds for col in axes(Z, 2)
        out[col], _ = VineCopulas._reduced_probability_integrand!(
            ws, plan, case.lower, case.upper, view(Z, :, col), order, slots,
            values,
        )
    end
    return out
end

function dense_ranks(values)
    permutation = sortperm(values; alg=Base.Sort.MergeSort)
    ranks = zeros(Float64, length(values))
    i = 1
    while i <= length(values)
        j = i
        while j < length(values) && values[permutation[j + 1]] == values[permutation[i]]
            j += 1
        end
        rank = (i + j) / 2
        for k in i:j
            ranks[permutation[k]] = rank
        end
        i = j + 1
    end
    return ranks
end

function spearman(x, y)
    rx, ry = dense_ranks(x), dense_ranks(y)
    return cor(rx, ry)
end

function kendall(x, y)
    concordance = 0
    total = 0
    for i in 1:(length(x) - 1), j in (i + 1):length(x)
        sx = sign(x[j] - x[i])
        sy = sign(y[j] - y[i])
        (sx == 0 || sy == 0) && continue
        concordance += sx * sy
        total += 1
    end
    return total == 0 ? NaN : concordance / total
end

function mixed_central_difference(f, z, h)
    s = length(z)
    total = 0.0
    point = similar(z)
    for mask in 0:(2^s - 1)
        coefficient = 1.0
        @inbounds for j in 1:s
            direction = ((mask >> (j - 1)) & 1) == 1 ? 1.0 : -1.0
            point[j] = z[j] + direction * h
            coefficient *= direction
        end
        total += coefficient * f(point)
    end
    return total / (2h)^s
end

function attempt_forwarddiff(case, plan, order, slots)
    z0 = fill(0.43, length(order))
    ws = VineCopulas._reduced_workspace(plan)
    values = fill(NaN, length(case.vine))
    f(x) = VineCopulas._reduced_probability_integrand!(
        ws, plan, case.lower, case.upper, x, order, slots, values,
    )[1]
    try
        VineCopulas.ForwardDiff.gradient(f, z0)
        return "available"
    catch err
        return "unavailable: " * replace(sprint(showerror, err), ',' => ';', '\n' => ' ')
    end
end

order_rows = NamedTuple[]
variance_rows = NamedTuple[]
iid_rows = NamedTuple[]
rqmc_rows = NamedTuple[]
roughness_rows = NamedTuple[]

for (case_index, case) in enumerate(cases())
    p = length(case.vine)
    plan = VineCopulas._compile_reduced_probability(case.vine)
    orders = VineCopulas._reduced_admissible_orders(plan)
    ref, ref_uncertainty = reference(case.id)
    @printf("%s: %d admissible orders\n", case.id, length(orders))

    for (order_index, order) in enumerate(orders)
        push!(order_rows, (;
            case_id=case.id,
            order_index,
            order=order_string(order),
            terminal_a=plan.terminal.a,
            terminal_b=plan.terminal.b,
            conditioning=order_string(plan.terminal.conditioning),
            interior=case.interior,
        ))
    end

    # Common random numbers across all orders for the high-precision iid
    # variance criterion.
    rng_variance = Random.MersenneTwister(BASE_SEED + 10_000case_index)
    Z_variance = rand(rng_variance, p - 2, IID_VARIANCE_N)
    variance_samples = Vector{Float64}(undef, IID_VARIANCE_N)
    for (order_index, order) in enumerate(orders)
        slots = VineCopulas._reduced_sampling_slots(plan, order)
        integrand_samples!(variance_samples, case, plan, order, slots, Z_variance)
        mean_g = mean(variance_samples)
        variance_g = var(variance_samples; corrected=true)
        variance_se = sqrt(
            max(mean((variance_samples .- mean_g).^4) -
                ((IID_VARIANCE_N - 3) / (IID_VARIANCE_N - 1)) * variance_g^2, 0.0) /
            IID_VARIANCE_N,
        )
        push!(variance_rows, (;
            case_id=case.id,
            order_index,
            order=order_string(order),
            N=IID_VARIANCE_N,
            mean_integrand=mean_g,
            reference=ref,
            reference_uncertainty=ref_uncertainty,
            variance=variance_g,
            variance_se,
            q1=case.upper[first(order)] - case.lower[first(order)],
        ))
    end

    # Independent iid replicates, again paired across orders.
    for rep in 1:IID_REPS
        seed = BASE_SEED + 100_000case_index + rep
        Z = rand(Random.MersenneTwister(seed), p - 2, IID_N)
        samples = Vector{Float64}(undef, IID_N)
        for (order_index, order) in enumerate(orders)
            slots = VineCopulas._reduced_sampling_slots(plan, order)
            integrand_samples!(samples, case, plan, order, slots, Z)
            estimate = mean(samples)
            push!(iid_rows, (;
                case_id=case.id,
                order_index,
                order=order_string(order),
                N=IID_N,
                rep,
                seed,
                estimate,
                error=estimate - ref,
                abs_error=abs(estimate - ref),
            ))
        end
    end

    # Scrambled Sobol replicates use the same seed for every order in a paired
    # cell, hence the same randomized point set.
    for power in RQMC_POWERS, rep in 1:RQMC_REPS
        N = 1 << power
        seed = BASE_SEED + 1_000_000case_index + 10_000power + rep
        for (order_index, order) in enumerate(orders)
            timed = @timed VineCopulas._rectprob_reduced_rqmc_with_order(
                case.vine, case.lower, case.upper, order;
                N, R=1, seed,
            )
            estimate = timed.value.estimate
            push!(rqmc_rows, (;
                case_id=case.id,
                order_index,
                order=order_string(order),
                N,
                rep,
                seed,
                estimate,
                error=estimate - ref,
                abs_error=abs(estimate - ref),
                runtime_s=timed.time,
            ))
        end
    end

    # The complete mixed derivative is only a partial RQMC diagnostic.  The
    # production scalar engine converts inputs to Float64, so ForwardDiff is
    # audited first.  Stable central differences are restricted to the single
    # strictly interior d=5 case and two step sizes.
    if case.interior
        rng_roughness = Random.MersenneTwister(BASE_SEED + 2_000_000case_index)
        margin = 0.04
        Z = margin .+ (1 - 2margin) .* rand(rng_roughness, p - 2, ROUGHNESS_N)
        for (order_index, order) in enumerate(orders)
            slots = VineCopulas._reduced_sampling_slots(plan, order)
            ad_status = attempt_forwarddiff(case, plan, order, slots)
            ws = VineCopulas._reduced_workspace(plan)
            values = fill(NaN, p)
            f(z) = VineCopulas._reduced_probability_integrand!(
                ws, plan, case.lower, case.upper, z, order, slots, values,
            )[1]
            for h in (0.002, 0.001)
                derivatives = [mixed_central_difference(f, view(Z, :, k), h)
                               for k in axes(Z, 2)]
                roughness = (1 - 2margin)^(p - 2) * mean(abs2, derivatives)
                push!(roughness_rows, (;
                    case_id=case.id,
                    order_index,
                    order=order_string(order),
                    derivative="top_order_only",
                    method="central_difference",
                    h,
                    interior_margin=margin,
                    N=ROUGHNESS_N,
                    roughness,
                    median_abs_derivative=median(abs.(derivatives)),
                    p90_abs_derivative=quantile(abs.(derivatives), 0.90),
                    max_abs_derivative=maximum(abs.(derivatives)),
                    ad_status,
                ))
            end
        end
    end
end

write_rows(joinpath(OUT, "admissible_orders.csv"), order_rows)
write_rows(joinpath(OUT, "iid_variance_criterion.csv"), variance_rows)
write_rows(joinpath(OUT, "iid_raw.csv"), iid_rows)
write_rows(joinpath(OUT, "rqmc_raw.csv"), rqmc_rows)
write_rows(joinpath(OUT, "top_derivative_diagnostic.csv"), roughness_rows)

summary_rows = NamedTuple[]
for vrow in variance_rows
    iid_group = filter(r -> r.case_id == vrow.case_id && r.order == vrow.order, iid_rows)
    rqmc_group = filter(r -> r.case_id == vrow.case_id && r.order == vrow.order &&
                             r.N == maximum(1 .<< RQMC_POWERS), rqmc_rows)
    iid_rmse = sqrt(mean(r.error^2 for r in iid_group))
    rqmc_rmse = sqrt(mean(r.error^2 for r in rqmc_group))
    rough = filter(r -> r.case_id == vrow.case_id && r.order == vrow.order &&
                        r.h == 0.001, roughness_rows)
    push!(summary_rows, (;
        case_id=vrow.case_id,
        order_index=vrow.order_index,
        order=vrow.order,
        V_pi=vrow.variance,
        iid_rmse,
        iid_N_rmse2_over_V=IID_N * iid_rmse^2 / vrow.variance,
        rqmc_N=maximum(1 .<< RQMC_POWERS),
        rqmc_rmse,
        median_rqmc_runtime_s=median(r.runtime_s for r in rqmc_group),
        top_derivative_roughness=isempty(rough) ? NaN : only(rough).roughness,
    ))
end
write_rows(joinpath(OUT, "ordering_summary.csv"), summary_rows)

association_rows = NamedTuple[]
for case in cases()
    rows = filter(r -> r.case_id == case.id, summary_rows)
    V = [r.V_pi for r in rows]
    iid_efficiency = [r.iid_rmse^2 for r in rows]
    rqmc_efficiency = [r.rqmc_rmse^2 for r in rows]
    rough = [r.top_derivative_roughness for r in rows]
    push!(association_rows, (;
        case_id=case.id,
        admissible_orders=length(rows),
        best_V_order=rows[argmin(V)].order,
        worst_V_order=rows[argmax(V)].order,
        V_worst_best_ratio=maximum(V) / minimum(V),
        best_iid_order=rows[argmin(iid_efficiency)].order,
        iid_worst_best_ratio=maximum(iid_efficiency) / minimum(iid_efficiency),
        spearman_V_iid=spearman(V, iid_efficiency),
        kendall_V_iid=kendall(V, iid_efficiency),
        best_rqmc_order=rows[argmin(rqmc_efficiency)].order,
        rqmc_worst_best_ratio=maximum(rqmc_efficiency) / minimum(rqmc_efficiency),
        spearman_V_rqmc=spearman(V, rqmc_efficiency),
        kendall_V_rqmc=kendall(V, rqmc_efficiency),
        spearman_roughness_rqmc=all(isfinite, rough) ? spearman(rough, rqmc_efficiency) : NaN,
        kendall_roughness_rqmc=all(isfinite, rough) ? kendall(rough, rqmc_efficiency) : NaN,
    ))
end
write_rows(joinpath(OUT, "ordering_associations.csv"), association_rows)

println("Ordering experiment complete: ", OUT)
