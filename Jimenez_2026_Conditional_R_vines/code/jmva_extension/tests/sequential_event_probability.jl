@testitem "Experimental sequential events – admissible orders and pilot" tags=[:Vine, :Experimental, :SequentialEvent, :Order] setup=[M] begin
    V = M.rvine5_general()
    plan = VineCopulas._compile_reduced_probability(V)
    orders = VineCopulas._reduced_admissible_orders(plan)
    @test length(orders) == 4
    @test all(length(order) == 3 for order in orders)
    @test all(sort(order) == sort(plan.terminal.conditioning) for order in orders)
    @test plan.default_order in orders

    lower = zeros(5)
    upper = [0.36, 0.61, 0.47, 0.73, 0.58]
    pilot = VineCopulas._select_reduced_order_pilot(
        V, lower, upper; M=512, seed=901, candidate_orders=orders,
    )
    @test pilot.order in orders
    @test length(pilot.variances) == length(orders)
    @test all(isfinite, pilot.variances)
    @test all(>=(0), pilot.variances)

    selected = VineCopulas._rectprob_reduced_pilot_rqmc(
        V, lower, upper; M=256, N=256, R=2,
        pilot_seed=902, final_seed=903,
    )
    @test selected.pilot.order in orders
    @test 0 <= selected.result.estimate <= 1
    @test_throws ArgumentError VineCopulas._rectprob_reduced_pilot_rqmc(
        V, lower, upper; M=32, N=32, R=1,
        pilot_seed=904, final_seed=904,
    )
end

@testitem "Experimental sequential events – rotated terminal CDF safety" tags=[:Vine, :Experimental, :SequentialEvent, :Survival] setup=[M] begin
    using Copulas
    using Distributions
    for flips in ((1,), (2,), (1, 2))
        rotated = SurvivalCopula(ClaytonCopula(2, 2.0), flips)
        @test VineCopulas._pair_cdf_is_safe(rotated)
        buffer = zeros(2)
        for (u, v) in ((0.2, 0.7), (0.8, 0.4))
            @test VineCopulas._pair_cdf(rotated, u, v, buffer) ≈ cdf(rotated, [u, v]) atol=1e-14
        end
    end
end

@testitem "Experimental sequential events – ordinal independence" tags=[:Vine, :Experimental, :SequentialEvent, :Ordinal] setup=[M] begin
    independent() = Copulas.IndependentCopula(2)
    V = DVineCopula(1:3, [(independent(), independent()), (independent(),)])
    patterns = ([1, 2, 3], [1, 3, 2], [2, 1, 3],
                [2, 3, 1], [3, 1, 2], [3, 2, 1])
    values = Float64[]
    dimensions = Int[]
    for pattern in patterns
        result = VineCopulas._ordinal_probability_rqmc(
            V, pattern; N=1 << 13, R=1, randomized=false,
        )
        push!(values, result.estimate)
        push!(dimensions, result.dimension)
        @test result.estimate ≈ 1 / 6 atol=8e-6
    end
    @test sum(values) ≈ 1 atol=8e-6
    # The terminal pair of this D-vine is (1,3): patterns with label 2 between
    # them are nonadjacent and use one dimension; the others use two.
    @test count(==(1), dimensions) == 2
    @test count(==(2), dimensions) == 4
end

@testitem "Experimental sequential events – all d=5 patterns" tags=[:Vine, :Experimental, :SequentialEvent, :Ordinal] setup=[M] begin
    function all_permutations(values)
        isempty(values) && return [Int[]]
        result = Vector{Vector{Int}}()
        for i in eachindex(values)
            remaining = [values[j] for j in eachindex(values) if j != i]
            for tail in all_permutations(remaining)
                push!(result, vcat(values[i], tail))
            end
        end
        return result
    end

    gaussian(rho) = GaussianCopula([1.0 rho; rho 1.0])
    V = DVineCopula(
        1:5,
        [Tuple(gaussian(0.45 / (1 + 0.12 * (k - 1))) for _ in 1:(5-k))
         for k in 1:4],
    )
    plan = VineCopulas._compile_reduced_probability(V)
    patterns = all_permutations(collect(1:5))
    @test length(patterns) == 120
    adjacent = count(patterns) do pattern
        ranks = VineCopulas._ordinal_ranks(pattern, 5)
        abs(ranks[plan.terminal.a] - ranks[plan.terminal.b]) == 1
    end
    @test adjacent == 2 * factorial(4)

    estimates = Float64[]
    for (index, pattern) in enumerate(patterns)
        result = VineCopulas._ordinal_probability_rqmc(
            V, pattern; N=1 << 10, R=2, seed=1200 + 10index,
        )
        @test result.estimate >= 0
        expected_dimension = begin
            ranks = VineCopulas._ordinal_ranks(pattern, 5)
            3 + Int(abs(ranks[plan.terminal.a] - ranks[plan.terminal.b]) == 1)
        end
        @test result.dimension == expected_dimension
        push!(estimates, result.estimate)
    end
    @test sum(estimates) ≈ 1 atol=4e-3
end

@testitem "Experimental sequential events – ordinal integrand invariant" tags=[:Vine, :Experimental, :SequentialEvent, :Invariant] setup=[M] begin
    using Random
    V = M.rvine5_general()
    plan = VineCopulas._compile_reduced_probability(V)
    order = first(VineCopulas._reduced_admissible_orders(plan))
    slots = VineCopulas._reduced_sampling_slots(plan, order)
    rng = M.stable_rng(991)
    for pattern in ([1, 3, 5, 2, 4], [1, 5, 2, 3, 4])
        ranks = VineCopulas._ordinal_ranks(pattern, 5)
        adjacent = abs(ranks[plan.terminal.a] - ranks[plan.terminal.b]) == 1
        for _ in 1:100
            ws = VineCopulas._reduced_workspace(plan)
            values = fill(NaN, 5)
            sampled = falses(5)
            G, q1 = VineCopulas._ordinal_probability_integrand!(
                ws, plan, pattern, ranks, rand(rng, 3 + Int(adjacent)),
                order, slots, values, sampled,
            )
            @test 0 <= G <= q1 <= 1
        end
    end
end

@testitem "Experimental sequential events – normalized entropy gradient" tags=[:Vine, :Experimental, :SequentialEvent, :Entropy] begin
    raw = [0.17, 0.29, 0.41, 0.23]
    result = VineCopulas._normalized_permutation_entropy_with_gradient(raw)
    @test sum(result.probabilities) ≈ 1.0 atol=1e-15
    @test result.entropy ≈ -sum(result.probabilities .* log.(result.probabilities)) /
        log(length(raw)) atol=1e-15

    h = 1e-6
    for i in eachindex(raw)
        plus = copy(raw)
        minus = copy(raw)
        plus[i] += h
        minus[i] -= h
        finite_difference = (
            VineCopulas._normalized_permutation_entropy_with_gradient(plus).entropy -
            VineCopulas._normalized_permutation_entropy_with_gradient(minus).entropy
        ) / (2h)
        @test result.gradient[i] ≈ finite_difference rtol=2e-7 atol=2e-10
    end

    # Entropy is invariant to common rescaling of the raw probabilities, so
    # its directional derivative in the radial direction must vanish.
    @test sum(result.gradient .* raw) ≈ 0.0 atol=2e-15
    @test_throws DomainError VineCopulas._normalized_permutation_entropy_with_gradient(
        [0.0, 0.5, 0.5],
    )
end
