using Dates
using Distributions
using Printf
using Serialization
using Statistics
using VineCopulas

include(joinpath(@__DIR__, "..", "..", "publication", "publication_common.jl"))

const DATA_FILE = joinpath(@__DIR__, "data", "10_Industry_Portfolios_daily_CSV.zip")
const OUT = get(ENV, "OUT", joinpath(@__DIR__, "results"))
const FIT_START = Date(get(ENV, "FIT_START", "2015-01-01"))
const FIT_END = Date(get(ENV, "FIT_END", "2024-12-31"))
const EMBEDDING = parse(Int, get(ENV, "EMBEDDING", "5"))
const DELAY = parse(Int, get(ENV, "DELAY", "1"))
const N = parse(Int, get(ENV, "N", string(1 << 12)))
const R = parse(Int, get(ENV, "R", "8"))
const BASE_SEED = parse(Int, get(ENV, "BASE_SEED", "840001"))
const INDUSTRY_FILTER = filter(!isempty, split(get(ENV, "INDUSTRIES", ""), ','))
const TRANSFORM_FILTER = filter(!isempty, split(get(ENV, "TRANSFORMS", "returns,absolute"), ','))

EMBEDDING == 5 || error("the frozen JMVA application design uses embedding dimension 5")
DELAY == 1 || error("the frozen JMVA application design uses delay 1")
ispow2(N) || error("N must be a power of two")
isfile(DATA_FILE) || error("missing fixed Kenneth French input: $DATA_FILE")
mkpath(OUT)
mkpath(joinpath(OUT, "models"))

function read_zip_text(path)
    return read(`unzip -p $path`, String)
end

function parse_industry_daily(text)
    lines = split(text, '\n')
    header_idx = findfirst(line -> startswith(strip(line), ",NoDur"), lines)
    header_idx === nothing && error("could not locate 10 Industry daily header")
    names = split(strip(lines[header_idx]), ',')[2:end]
    dates = Date[]
    rows = Vector{Float64}[]
    for line in Iterators.drop(lines, header_idx)
        fields = split(strip(line), ',')
        length(fields) == 11 || continue
        occursin(r"^\d{8}$", fields[1]) || continue
        values = try
            parse.(Float64, fields[2:end]) ./ 100
        catch
            continue
        end
        any(x -> x <= -0.99 || !isfinite(x), values) && continue
        push!(dates, Date(fields[1], dateformat"yyyymmdd"))
        push!(rows, values)
    end
    X = reduce(vcat, permutedims.(rows))
    return names, dates, X
end

function midrank_uniform(x)
    n = length(x)
    order = sortperm(x; alg=Base.Sort.MergeSort)
    ranks = zeros(Float64, n)
    i = 1
    while i <= n
        j = i
        while j < n && x[order[j + 1]] == x[order[i]]
            j += 1
        end
        rank = (i + j) / 2
        @inbounds for k in i:j
            ranks[order[k]] = rank
        end
        i = j + 1
    end
    return ranks ./ (n + 1)
end

function lag_matrix(u, m, delay)
    n = length(u) - (m - 1) * delay
    n > 0 || error("series is too short for requested embedding")
    X = Matrix{Float64}(undef, m, n)
    @inbounds for col in 1:n, row in 1:m
        X[row, col] = u[col + (row - 1) * delay]
    end
    return X
end

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

function ordinal_pattern(x)
    return sortperm(collect(eachindex(x)); by=i -> (x[i], i), alg=Base.Sort.MergeSort)
end

function empirical_patterns(x, patterns, m, delay)
    keys = Dict(join(pattern, '-') => i for (i, pattern) in enumerate(patterns))
    counts = zeros(Int, length(patterns))
    tied = 0
    windows = length(x) - (m - 1) * delay
    for start in 1:windows
        window = [x[start + (j - 1) * delay] for j in 1:m]
        if length(unique(window)) < m
            tied += 1
            continue
        end
        key = join(ordinal_pattern(window), '-')
        counts[keys[key]] += 1
    end
    complete = windows - tied
    complete > 0 || error("no tie-free ordinal windows")
    return counts ./ complete, tied, windows, complete
end

function js_divergence(p, q)
    midpoint = (p .+ q) ./ 2
    term(x, y) = iszero(x) ? 0.0 : x * log(x / y)
    return 0.5sum(term(p[i], midpoint[i]) for i in eachindex(p)) +
           0.5sum(term(q[i], midpoint[i]) for i in eachindex(q))
end

function family_description(vine)
    return join((VineCopulas._short_family_name(edge.copula)
                 for edge in VineCopulas.vine_edges(vine)), ';')
end

names, dates, Xall = parse_industry_daily(read_zip_text(DATA_FILE))
keep = findall(date -> FIT_START <= date <= FIT_END, dates)
X = Xall[keep, :]
dates_kept = dates[keep]
patterns = all_permutations(collect(1:EMBEDDING))

raw_rows = NamedTuple[]
model_rows = NamedTuple[]
summary_rows = NamedTuple[]

for (industry_index, industry) in enumerate(names)
    !isempty(INDUSTRY_FILTER) && industry ∉ INDUSTRY_FILTER && continue
    base_returns = collect(view(X, :, industry_index))
    for (transform_index, transform) in enumerate(("returns", "absolute"))
        transform ∉ TRANSFORM_FILTER && continue
        series = transform == "returns" ? base_returns : abs.(base_returns)
        empirical, tied_windows, windows, complete_windows = empirical_patterns(
            series, patterns, EMBEDDING, DELAY,
        )

        # A single common empirical margin is applied before lagging. This
        # preserves the ordinal comparisons and reflects the common-margin
        # temporal-copula formulation.
        U = lag_matrix(midrank_uniform(series), EMBEDDING, DELAY)
        fit_timed = @timed fit(
            RVineCopula, U;
            family_set=(GaussianCopula, ClaytonCopula, GumbelCopula,
                        FrankCopula, JoeCopula),
            include_independence=true,
            allow_rotations=true,
            preselect=true,
            selection_criterion=:bic,
            tree_criterion=:tau,
            trace=false,
        )
        vine = fit_timed.value
        model_file = joinpath(OUT, "models", "$(industry)_$(transform).jls")
        serialize(model_file, vine)
        plan = VineCopulas._compile_reduced_probability(vine)
        push!(model_rows, (;
            industry,
            transform,
            observations=length(series),
            lag_vectors=size(U, 2),
            first_date=first(dates_kept),
            last_date=last(dates_kept),
            tied_windows,
            fitting_time_s=fit_timed.time,
            selection_criterion="bic",
            tree_criterion="absolute_kendall_tau",
            candidate_families="Gaussian;Clayton;Gumbel;Frank;Joe;Independence",
            families=family_description(vine),
            terminal_a=plan.terminal.a,
            terminal_b=plan.terminal.b,
            admissible_orders=length(VineCopulas._reduced_admissible_orders(plan)),
            model_file=relpath(model_file, OUT),
        ))

        probabilities = zeros(Float64, length(patterns))
        numerical_se = zeros(Float64, length(patterns))
        for (pattern_index, pattern) in enumerate(patterns)
            seed = BASE_SEED + 1_000_000industry_index +
                   100_000transform_index + 100pattern_index
            timed = @timed VineCopulas._ordinal_probability_rqmc(
                vine, pattern; N, R, seed,
            )
            result = timed.value
            probabilities[pattern_index] = result.estimate
            numerical_se[pattern_index] = result.stderr
            ranks = VineCopulas._ordinal_ranks(pattern, EMBEDDING)
            adjacent = abs(ranks[plan.terminal.a] - ranks[plan.terminal.b]) == 1
            push!(raw_rows, (;
                industry,
                transform,
                pattern_index,
                pattern=join(pattern, '-'),
                terminal_adjacent=adjacent,
                integration_dimension=result.dimension,
                N,
                R,
                seed,
                model_probability=result.estimate,
                numerical_se=result.stderr,
                empirical_frequency=empirical[pattern_index],
                runtime_s=timed.time,
            ))
        end

        mass = sum(probabilities)
        mass > 0 || error("zero ordinal probability mass for $industry $transform")
        entropy_result = VineCopulas._normalized_permutation_entropy_with_gradient(
            probabilities,
        )
        normalized = entropy_result.probabilities
        entropy_model = entropy_result.entropy
        entropy_empirical = -sum(x -> iszero(x) ? 0.0 : x * log(x), empirical) /
            log(length(patterns))
        entropy_numerical_se = sqrt(sum(
            (entropy_result.gradient[i] * numerical_se[i])^2
            for i in eachindex(normalized)
        ))
        l1 = sum(abs, normalized .- empirical)
        push!(summary_rows, (;
            industry,
            transform,
            observations=length(series),
            lag_vectors=windows,
            empirical_tie_free_windows=complete_windows,
            tied_windows,
            tied_fraction=tied_windows / windows,
            model_mass=mass,
            mass_error=mass - 1,
            mass_numerical_se=sqrt(sum(abs2, numerical_se)),
            minimum_probability=minimum(probabilities),
            maximum_probability=maximum(probabilities),
            permutation_entropy_model=entropy_model,
            permutation_entropy_numerical_se=entropy_numerical_se,
            permutation_entropy_empirical=entropy_empirical,
            entropy_difference=entropy_model - entropy_empirical,
            l1_distance=l1,
            total_variation=0.5l1,
            jensen_shannon=js_divergence(normalized, empirical),
            fitting_time_s=fit_timed.time,
            probability_time_s=sum(r.runtime_s for r in raw_rows
                                   if r.industry == industry && r.transform == transform),
        ))
        @printf("%-6s %-8s H_model=%.6f H_emp=%.6f mass=%.8f\n",
                industry, transform, entropy_model, entropy_empirical, mass)
    end
end

write_rows(joinpath(OUT, "model_metadata.csv"), model_rows)
write_rows(joinpath(OUT, "ordinal_probabilities.csv"), raw_rows)
write_rows(joinpath(OUT, "permutation_entropy_summary.csv"), summary_rows)

open(joinpath(OUT, "DATA_PROVENANCE.md"), "w") do io
    println(io, "# Data provenance")
    println(io)
    println(io, "Source: Kenneth R. French Data Library, 10 Industry Portfolios, daily value-weighted returns.")
    println(io)
    println(io, "- Source URL: https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/10_Industry_Portfolios_daily_CSV.zip")
    println(io, "- Fixed local input: `data/10_Industry_Portfolios_daily_CSV.zip`")
    println(io, "- Analysis dates: $FIT_START through $FIT_END")
    println(io, "- Retained observations: $(length(keep))")
    println(io, "- First/last retained dates: $(first(dates_kept)) / $(last(dates_kept))")
    println(io, "- Embedding: m=$EMBEDDING, delay=$DELAY, overlapping windows")
    println(io, "- Common marginal transform: mid-ranks over each complete series divided by n+1, followed by lagging")
    println(io, "- Rounded ties and model fit: mid-ranks are a practical latent-continuous rank approximation, not a claim of exact observed continuity")
    println(io, "- Model selection: BIC is used as a working/composite pair-family score on overlapping lag-vector rows; no iid likelihood-based selection inference is claimed")
    println(io, "- Empirical ties: tied ordinal windows are excluded and reported in the summary; no iid inference is attached to overlapping frequencies")
end

println("Permutation-entropy application complete: ", OUT)
