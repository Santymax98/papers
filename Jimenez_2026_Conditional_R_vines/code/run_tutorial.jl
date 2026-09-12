using Printf
using Copulas
using VineCopulas

"""Return the lower/upper rectangle probability of the columns of `U`."""
function rectangle_indicator(U, lower, upper)
    count = 0
    @inbounds for j in axes(U, 2)
        inside = true
        for i in axes(U, 1)
            if !(lower[i] < U[i, j] <= upper[i])
                inside = false
                break
            end
        end
        count += inside
    end
    return count / size(U, 2)
end

gaussian_pair(rho) = GaussianCopula([1.0 rho; rho 1.0])

# A fixed, simplified four-dimensional D-vine used only as a fast smoke test.
vine = DVineCopula(
    [1, 2, 3, 4],
    [
        (gaussian_pair(0.35), gaussian_pair(0.55), gaussian_pair(-0.25)),
        (gaussian_pair(0.40), gaussian_pair(0.30)),
        (gaussian_pair(0.45),),
    ],
)

lower = zeros(4)
upper = [0.40, 0.55, 0.45, 0.60]
N = 1 << 10
seed = UInt64(20260906)

# Frozen full-dimensional baseline: Owen-scrambled Sobol points, inverse
# Rosenblatt transformation, and the rectangle indicator.
Z = VineCopulas._reduced_owen_points(length(vine), N, seed)
U = inverse_rosenblatt(vine, Z)
indicator_estimate = rectangle_indicator(U, lower, upper)

# Research implementation: truncated conditional transport with terminal-pair
# integration.  This private API is intentionally exercised directly here.
reduced = VineCopulas._rectprob_reduced_rqmc(
    vine, lower, upper; N=N, R=1, seed=Int(seed), randomized=true,
)

println("Conditional Dimension Reduction for Vine Copula Probabilities")
println("Smoke-test model: simplified Gaussian D-vine, d = 4")
println("Event: U_i <= upper_i, i = 1,...,4")
@printf("Budget N = %d; Owen-scramble seed = %d\n", N, seed)
@printf("Full-dimensional indicator estimate : %.12g\n", indicator_estimate)
@printf("Reduced conditional estimate        : %.12g\n", reduced.estimate)
println("Reduced RQMC replicate SE           : not defined for one replicate")
@printf("Absolute difference                 : %.4g\n", abs(indicator_estimate - reduced.estimate))

# JMVA extension: all six ordinal patterns of a three-dimensional independent
# temporal vine.  Their probabilities must sum to one and each equals 1/6.
independent_pair() = Copulas.IndependentCopula(2)
ordinal_vine = DVineCopula(
    1:3,
    [(independent_pair(), independent_pair()), (independent_pair(),)],
)
patterns = ([1, 2, 3], [1, 3, 2], [2, 1, 3],
            [2, 3, 1], [3, 1, 2], [3, 2, 1])
ordinal_probabilities = [
    VineCopulas._ordinal_probability_rqmc(
        ordinal_vine, pattern; N=1 << 11, R=1,
        seed=20261000 + i, randomized=false,
    ).estimate
    for (i, pattern) in enumerate(patterns)
]
@printf("Ordinal-pattern mass (d = 3)       : %.12g\n", sum(ordinal_probabilities))
@printf("Maximum error from 1/6             : %.4g\n",
        maximum(abs.(ordinal_probabilities .- 1 / 6)))
@assert isapprox(sum(ordinal_probabilities), 1.0; atol=5e-5)
println("Smoke test completed successfully.")
