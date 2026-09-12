using Copulas
using Dates
using Distributions
using Printf
using Random
using Statistics
using TOML
using VineCopulas

const PUBLICATION_DIR = @__DIR__
const REPO_ROOT = normpath(joinpath(PUBLICATION_DIR, "..", ".."))

_env_string(name, default) = get(ENV, name, default)
_env_int(name, default) = parse(Int, get(ENV, name, string(default)))
_env_float(name, default) = parse(Float64, get(ENV, name, string(default)))
_env_bool(name, default=false) = lowercase(get(ENV, name, default ? "true" : "false")) in ("1", "true", "t", "yes", "y")
_env_ints(name, default) = parse.(Int, split(get(ENV, name, default), ','))

function publication_config()
    mode = lowercase(_env_string("MODE", "pilot"))
    mode in ("pilot", "final") || error("MODE must be pilot or final")
    default_reps = mode == "final" ? 30 : 4
    default_ref_r = mode == "final" ? 32 : 8
    default_ref_powers = mode == "final" ? "16,18,20" : "14,16"
    default_bootstrap = mode == "final" ? 2000 : 300
    return (;
        mode,
        case_set=lowercase(_env_string("CASE_SET", "core")),
        dims=_env_ints("DIMS", "4,5,10,20"),
        powers=_env_ints("POWERS", "8,10,12,14,16"),
        reps=_env_int("REPS", default_reps),
        base_seed=_env_int("BASE_SEED", 271828),
        reference_seed=_env_int("REFERENCE_BASE_SEED", 8_675_309),
        reference_r=_env_int("REFERENCE_R", default_ref_r),
        reference_powers=_env_ints("REFERENCE_POWERS", default_ref_powers),
        reference_target_se=_env_float("REFERENCE_TARGET_SE", mode == "final" ? 1e-7 : 5e-6),
        reference_split_z_max=_env_float("REFERENCE_SPLIT_Z_MAX", 4.0),
        reference_maxevals=_env_int("REFERENCE_MAXEVALS", mode == "final" ? 5_000_000 : 1_000_000),
        reference_hcubature_rtol=_env_float("REFERENCE_HCUBATURE_RTOL", mode == "final" ? 2e-10 : 2e-9),
        reference_hcubature_atol=_env_float("REFERENCE_HCUBATURE_ATOL", mode == "final" ? 2e-12 : 2e-11),
        max_reference_ratio=_env_float("MAX_REFERENCE_RMSE_RATIO", 0.20),
        bootstrap_b=_env_int("BOOTSTRAP_B", default_bootstrap),
        bootstrap_seed=_env_int("BOOTSTRAP_SEED", 314159),
        force_references=_env_bool("FORCE_REFERENCES", false),
        resume=_env_bool("RESUME", true),
        out=_env_string("OUT", joinpath("benchmarks", "results", "publication", Dates.format(now(), "yyyymmdd_HHMMSS"))),
    )
end

# -----------------------------
# Publication cases
# -----------------------------

gaussian_pair(rho) = GaussianCopula([1.0 rho; rho 1.0])

function publication_homogeneous_dvine(d::Int, rho::Real)
    return DVineCopula(
        collect(1:d),
        [Tuple(gaussian_pair(rho / (1 + 0.12 * (k - 1))) for _ in 1:(d-k)) for k in 1:(d-1)],
    )
end

function publication_mixed_dvine(d::Int)
    pool = PairCopula[
        gaussian_pair(0.45), ClaytonCopula(2, 1.4), GumbelCopula(2, 1.3),
        FrankCopula(2, 2.5), JoeCopula(2, 1.35),
    ]
    levels = [Tuple(pool[mod1(3k + i, length(pool))] for i in 1:(d-k)) for k in 1:(d-1)]
    # Keep the terminal pair on a deterministic analytical family.
    levels[end] = (FrankCopula(2, 2.5),)
    return DVineCopula(collect(1:d), levels)
end

function publication_genuine_rvine5()
    ord = [1, 3, 2, 4, 5]
    S = ([2, 2, 4, 5], [3, 4, 5], [4, 5], [5])
    E = (
        (gaussian_pair(0.35), ClaytonCopula(2, 1.4), FrankCopula(2, 2.0), GumbelCopula(2, 1.25)),
        (FrankCopula(2, 1.7), gaussian_pair(-0.25), ClaytonCopula(2, 1.1)),
        (GumbelCopula(2, 1.3), gaussian_pair(0.20)),
        (ClaytonCopula(2, 0.8),),
    )
    return RVineCopula(ord, S, E)
end

function publication_genuine_rvine7()
    ord = [4, 3, 7, 1, 2, 5, 6]
    S = ([5, 2, 6, 6, 6, 6], [6, 6, 1, 2, 5], [2, 5, 2, 5],
         [1, 1, 5], [3, 7], [7])
    pool = PairCopula[
        gaussian_pair(0.35), ClaytonCopula(2, 1.25), FrankCopula(2, 2.0),
        GumbelCopula(2, 1.2), JoeCopula(2, 1.3),
    ]
    E = [Tuple(pool[mod1(2k + i, length(pool))] for i in 1:(7-k)) for k in 1:6]
    E[end] = (FrankCopula(2, 2.2),)
    return RVineCopula(ord, S, E)
end

function publication_cases(cfg=publication_config())
    cases = NamedTuple[]
    dims = sort(unique(cfg.dims))

    for d in dims
        V = publication_homogeneous_dvine(d, 0.45)
        push!(cases, (;
            id="gaussian_moderate_d$d", design="gaussian_moderate", model="D",
            point="center", vine=V, lower=zeros(d), upper=fill(0.55, d),
            paper_group="main",
        ))
        push!(cases, (;
            id="gaussian_tail_d$d", design="gaussian_moderate", model="D",
            point="lower_asymmetric", vine=V, lower=zeros(d),
            upper=collect(range(0.12, 0.72; length=d)), paper_group="main",
        ))
        push!(cases, (;
            id="gaussian_rectangle_d$d", design="gaussian_moderate", model="D",
            point="rectangle", vine=V, lower=fill(0.08, d),
            upper=collect(range(0.48, 0.88; length=d)), paper_group="main",
        ))
        push!(cases, (;
            id="gaussian_strong_d$d", design="gaussian_strong", model="D",
            point="upper", vine=publication_homogeneous_dvine(d, 0.82),
            lower=zeros(d), upper=fill(0.72, d), paper_group="main",
        ))
    end

    for d in intersect(dims, [5, 10, 20])
        push!(cases, (;
            id="mixed_d$d", design="mixed", model="D", point="asymmetric",
            vine=publication_mixed_dvine(d), lower=zeros(d),
            upper=collect(range(0.25, 0.83; length=d)), paper_group="main",
        ))
    end

    if 5 in dims
        push!(cases, (;
            id="genuine_rvine_d5", design="genuine_rvine", model="R",
            point="asymmetric", vine=publication_genuine_rvine5(), lower=zeros(5),
            upper=[0.36, 0.61, 0.47, 0.73, 0.58], paper_group="main",
        ))
    end

    # The 7D genuine R-vine is useful for supplementary structural/order checks.
    if cfg.case_set == "extended" && 7 in dims
        push!(cases, (;
            id="genuine_rvine_d7", design="genuine_rvine", model="R",
            point="rectangle", vine=publication_genuine_rvine7(), lower=fill(0.04, 7),
            upper=[0.42, 0.68, 0.51, 0.76, 0.59, 0.81, 0.63], paper_group="supplement",
        ))
    end

    if cfg.case_set == "compact"
        keep = Set(["gaussian_moderate", "gaussian_strong"])
        cases = filter(c -> c.design in keep && c.point != "rectangle", cases)
    elseif !(cfg.case_set in ("core", "extended"))
        error("CASE_SET must be compact, core, or extended")
    end
    return cases
end

# -----------------------------
# CSV / metadata utilities
# -----------------------------

_csv_string(x::AbstractFloat) = isfinite(x) ? repr(Float64(x)) : string(x)
_csv_string(x) = string(x)

# Minimal RFC 4180-compatible CSV support.  Publication metadata can legitimately
# contain commas (e.g. Julia parametric type names), quotes, or newlines, so the
# benchmark must preserve those fields rather than silently sanitizing them.
function _csv_escape(x::AbstractString)
    if occursin(',', x) || occursin('"', x) || occursin('\n', x) || occursin('\r', x)
        return "\"" * replace(x, "\"" => "\"\"") * "\""
    end
    return x
end

function _csv_parse_line(line::AbstractString)
    fields = String[]
    buf = IOBuffer()
    in_quotes = false
    i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if in_quotes
            if c == '"'
                ni = nextind(line, i)
                if ni <= lastindex(line) && line[ni] == '"'
                    write(buf, '"')
                    i = ni
                else
                    in_quotes = false
                end
            else
                write(buf, c)
            end
        elseif c == '"'
            in_quotes = true
        elseif c == ','
            push!(fields, String(take!(buf)))
        else
            write(buf, c)
        end
        i = nextind(line, i)
    end
    in_quotes && error("unterminated quoted CSV field")
    push!(fields, String(take!(buf)))
    return fields
end

function write_rows(path::AbstractString, rows)
    isempty(rows) && error("cannot write empty CSV: $path")
    mkpath(dirname(path))
    open(path, "w") do io
        names = propertynames(first(rows))
        println(io, join((_csv_escape(string(n)) for n in names), ','))
        for row in rows
            fields = (_csv_escape(_csv_string(getproperty(row, n))) for n in names)
            println(io, join(fields, ','))
        end
    end
end

function append_row(path::AbstractString, row)
    mkpath(dirname(path))
    newfile = !isfile(path) || filesize(path) == 0
    open(path, "a") do io
        names = propertynames(row)
        newfile && println(io, join((_csv_escape(string(n)) for n in names), ','))
        fields = (_csv_escape(_csv_string(getproperty(row, n))) for n in names)
        println(io, join(fields, ','))
        flush(io)
    end
end

function read_table(path::AbstractString)
    isfile(path) || return Dict{String,String}[]
    lines = readlines(path)
    isempty(lines) && return Dict{String,String}[]
    header = _csv_parse_line(first(lines))
    rows = Dict{String,String}[]
    for line in Iterators.drop(lines, 1)
        isempty(strip(line)) && continue
        fields = _csv_parse_line(line)
        length(fields) == length(header) || error("malformed CSV row in $path")
        push!(rows, Dict(header .=> fields))
    end
    return rows
end

_parse_float(row, key) = parse(Float64, row[key])
_parse_int(row, key) = parse(Int, row[key])
_parse_bool(row, key) = lowercase(row[key]) == "true"
_vector_field(x) = join((repr(Float64(v)) for v in x), '|')

function _shell_output(args::Vector{String}; dir=REPO_ROOT)
    try
        return cd(dir) do
            strip(read(Cmd(args), String))
        end
    catch err
        @warn "Unable to capture shell metadata" command=join(args, " ") directory=dir exception=(err, catch_backtrace())
        return "unavailable"
    end
end

function capture_environment!(cfg, cases)
    out = cfg.out
    mkpath(out)
    envdir = joinpath(out, "environment")
    mkpath(envdir)
    for file in ("Project.toml", "Manifest.toml")
        src = joinpath(REPO_ROOT, "benchmarks", file)
        isfile(src) && cp(src, joinpath(envdir, file); force=true)
    end
    package_project = joinpath(REPO_ROOT, "Project.toml")
    isfile(package_project) && cp(package_project, joinpath(envdir, "package_Project.toml"); force=true)

    meta = Dict{String,Any}(
        "created_at" => string(now()),
        "julia_version" => string(VERSION),
        "machine" => Sys.MACHINE,
        "cpu_name" => Sys.CPU_NAME,
        "threads" => Threads.nthreads(),
        "word_size" => Sys.WORD_SIZE,
        "git_commit" => _shell_output(["git", "rev-parse", "HEAD"]),
        "git_branch" => _shell_output(["git", "branch", "--show-current"]),
        "git_status_porcelain" => _shell_output(["git", "status", "--porcelain"]),
        "mode" => cfg.mode,
        "case_set" => cfg.case_set,
        "dimensions" => cfg.dims,
        "powers" => cfg.powers,
        "reps" => cfg.reps,
        "base_seed" => cfg.base_seed,
        "reference_seed" => cfg.reference_seed,
        "reference_r" => cfg.reference_r,
        "reference_powers" => cfg.reference_powers,
        "reference_target_se" => cfg.reference_target_se,
        "reference_split_z_max" => cfg.reference_split_z_max,
        "max_reference_rmse_ratio" => cfg.max_reference_ratio,
        "bootstrap_b" => cfg.bootstrap_b,
        "n_cases" => length(cases),
        "rqmc_rule" => "Sobol + Owen nested uniform scrambling; independent deterministic axis seeds",
        "current_method" => "d-dimensional inverse Rosenblatt followed by rectangle indicator",
        "reduced_method" => "(d-2)-dimensional sequentially truncated inverse Rosenblatt + terminal pair probability",
    )
    open(joinpath(out, "run_metadata.toml"), "w") do io
        TOML.print(io, meta)
    end

    manifest_rows = NamedTuple[]
    for case in cases
        plan = VineCopulas._compile_reduced_probability(case.vine)
        push!(manifest_rows, (;
            case_id=case.id,
            design=case.design,
            model=case.model,
            point=case.point,
            paper_group=case.paper_group,
            d=length(case.vine),
            conditioning_dimension=length(case.vine)-2,
            lower=_vector_field(case.lower),
            upper=_vector_field(case.upper),
            is_rectangle=any(x -> !iszero(x), case.lower),
            terminal_pair=string(typeof(plan.terminal_copula)),
            default_order=join(plan.default_order, '-'),
        ))
    end
    write_rows(joinpath(out, "case_manifest.csv"), manifest_rows)
    return nothing
end

# -----------------------------
# Estimators
# -----------------------------

function rectangle_indicator(U, lower, upper)
    count = 0
    @inbounds for col in axes(U, 2)
        inside = true
        for row in axes(U, 1)
            if !(lower[row] < U[row, col] <= upper[row])
                inside = false
                break
            end
        end
        count += inside
    end
    return count / size(U, 2)
end

function current_estimate(vine, lower, upper, N::Int, seed::Int)
    Z = VineCopulas._reduced_owen_points(length(vine), N, UInt64(seed))
    U = inverse_rosenblatt(vine, Z)
    return rectangle_indicator(U, lower, upper)
end

mutable struct PublicationReducedRunner{P,W}
    plan::P
    order::Vector{Int}
    slots::Vector{Int}
    workspace::W
    conditioning_values::Vector{Float64}
    lower::Vector{Float64}
    upper::Vector{Float64}
end

function prepare_reduced_runner(case; order_strategy::Symbol=:default)
    p = length(case.vine)
    lo, hi = VineCopulas._check_reduced_bounds(p, case.lower, case.upper)
    plan = VineCopulas._compile_reduced_probability(case.vine)
    order = VineCopulas._reduced_sampling_order(plan, hi .- lo, order_strategy)
    slots = VineCopulas._reduced_sampling_slots(plan, order)
    workspace = VineCopulas._reduced_workspace(plan)
    values = fill(NaN, p)
    return PublicationReducedRunner(plan, order, slots, workspace, values, lo, hi)
end

function reset_counts!(runner::PublicationReducedRunner)
    runner.workspace.counts.hfunc = 0
    runner.workspace.counts.hinv = 0
    runner.workspace.counts.pair_cdf = 0
    return nothing
end

function reduced_compiled_estimate!(runner::PublicationReducedRunner, N::Int, seed::Int)
    reset_counts!(runner)
    p = runner.plan.p
    if any(i -> runner.lower[i] == runner.upper[i], 1:p)
        return (; estimate=0.0, hfunc_calls=0, hinv_calls=0, pair_cdf_calls=0)
    end
    if all(iszero, runner.lower) && all(isone, runner.upper)
        return (; estimate=1.0, hfunc_calls=0, hinv_calls=0, pair_cdf_calls=0)
    end
    Z = VineCopulas._reduced_owen_points(p - 2, N, UInt64(seed))
    total = 0.0
    @inbounds for col in axes(Z, 2)
        G, _ = VineCopulas._reduced_probability_integrand!(
            runner.workspace,
            runner.plan,
            runner.lower,
            runner.upper,
            view(Z, :, col),
            runner.order,
            runner.slots,
            runner.conditioning_values,
        )
        total += G
    end
    return (;
        estimate=total / N,
        hfunc_calls=runner.workspace.counts.hfunc,
        hinv_calls=runner.workspace.counts.hinv,
        pair_cdf_calls=runner.workspace.counts.pair_cdf,
    )
end

_gc_allocs(timed) = Base.gc_alloc_count(timed.gcstats)

# -----------------------------
# Reference helpers
# -----------------------------

reference_dir(cfg) = joinpath(cfg.out, "references")
reference_toml(cfg, case_id) = joinpath(reference_dir(cfg), case_id * ".toml")

function reference_integrand(case)
    p = length(case.vine)
    runner = prepare_reduced_runner(case)
    f = z -> begin
        G, _ = VineCopulas._reduced_probability_integrand!(
            runner.workspace,
            runner.plan,
            runner.lower,
            runner.upper,
            z,
            runner.order,
            runner.slots,
            runner.conditioning_values,
        )
        G
    end
    return f, p - 2
end

function split_half_diagnostic(replicates::AbstractVector{<:Real})
    R = length(replicates)
    R >= 4 || return (NaN, NaN, NaN, NaN, NaN)
    h = R ÷ 2
    a = Float64.(replicates[1:h])
    b = Float64.(replicates[(h+1):end])
    ma, mb = mean(a), mean(b)
    sea = length(a) > 1 ? std(a) / sqrt(length(a)) : NaN
    seb = length(b) > 1 ? std(b) / sqrt(length(b)) : NaN
    denom = sqrt(sea^2 + seb^2)
    z = denom > 0 ? abs(ma - mb) / denom : 0.0
    return ma, mb, sea, seb, z
end

function load_reference(cfg, case_id::AbstractString)
    path = reference_toml(cfg, case_id)
    isfile(path) || error("missing reference for $case_id; run generate_references.jl first")
    return TOML.parsefile(path)
end

function aggregate_reference_files!(cfg, cases)
    rows = NamedTuple[]
    for case in cases
        ref = load_reference(cfg, case.id)
        push!(rows, (;
            case_id=case.id,
            d=length(case.vine),
            method=ref["method"],
            estimate=ref["estimate"],
            uncertainty=ref["uncertainty"],
            uncertainty_kind=ref["uncertainty_kind"],
            N=ref["N"],
            R=ref["R"],
            seed=ref["seed"],
            target_met=ref["target_met"],
            split_half_z=ref["split_half_z"],
        ))
    end
    write_rows(joinpath(cfg.out, "references.csv"), rows)
    return rows
end
