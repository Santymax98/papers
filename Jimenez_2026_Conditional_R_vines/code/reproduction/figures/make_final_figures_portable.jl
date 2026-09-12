using CSV
using DataFrames
using CairoMakie
using Statistics

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const FIGDIR = joinpath(@__DIR__, "outputs")
mkpath(FIGDIR)

const CHECKPOINT = joinpath(ROOT, "internal_benchmark", "frozen_results")
const INTEGRATION_SUMMARY = joinpath(ROOT, "reproduction", "integration_rules", "summaries", "integration_rule_summary.csv")
const HIGHDIM_SUMMARY = joinpath(ROOT, "reproduction", "highdim", "summaries", "highdim_summary.csv")
const FINANCIAL_SUMMARY = joinpath(ROOT, "reproduction", "financial", "summaries", "05_kenneth_french_summary.csv")

set_theme!(
    Theme(
        fontsize = 11,
        Axis = (
            xgridvisible = true,
            ygridvisible = true,
            xgridcolor = (:gray, 0.18),
            ygridcolor = (:gray, 0.18),
            xlabelsize = 11,
            ylabelsize = 11,
            titlesize = 11,
        ),
        Legend = (framevisible = false, labelsize = 10),
    ),
)

const METHOD_LABEL = Dict(
    "current_indicator" => "Indicator",
    "reduced_conditional" => "Reduced",
)

const METHOD_COLOR = Dict(
    "current_indicator" => :gray35,
    "reduced_conditional" => :dodgerblue3,
)

const METHOD_MARKER = Dict(
    "current_indicator" => :circle,
    "reduced_conditional" => :diamond,
)

function add_lines!(ax, df, xcol, ycol; label_suffix = "")
    for method in ["current_indicator", "reduced_conditional"]
        m = sort(df[df.method .== method, :], xcol)
        isempty(m) && continue
        lines!(ax, Float64.(m[!, xcol]), Float64.(m[!, ycol]);
            color = METHOD_COLOR[method], linewidth = 2,
            label = METHOD_LABEL[method] * label_suffix)
        scatter!(ax, Float64.(m[!, xcol]), Float64.(m[!, ycol]);
            color = METHOD_COLOR[method], marker = METHOD_MARKER[method],
            markersize = 8)
    end
end

function fig1_internal_rmse_runtime()
    df = CSV.read(joinpath(CHECKPOINT, "summary_by_case.csv"), DataFrame)
    cases = [
        "gaussian_moderate_d20",
        "gaussian_strong_d20",
        "gaussian_rectangle_d10",
        "mixed_d20",
        "genuine_rvine_d5",
    ]
    titles = Dict(
        "gaussian_moderate_d20" => "Gaussian moderate, d=20",
        "gaussian_strong_d20" => "Gaussian strong, d=20",
        "gaussian_rectangle_d10" => "Interior rectangle, d=10",
        "mixed_d20" => "Mixed families, d=20",
        "genuine_rvine_d5" => "Genuine R-vine, d=5",
    )
    fig = Figure(size = (720, 520))
    positions = [(1,1), (1,2), (1,3), (2,1), (2,2)]
    for (case, pos) in zip(cases, positions)
        sub = df[(df.case_id .== case) .&& in.(df.method, Ref(["current_indicator", "reduced_conditional"])), :]
        ax = Axis(fig[pos...], title = titles[case], xscale = log10, yscale = log10,
            xlabel = "median runtime (s)", ylabel = "RMSE")
        add_lines!(ax, sub, :median_runtime_s, :rmse)
        pos == positions[1] && axislegend(ax, position = :lb, unique = true)
    end
    save(joinpath(FIGDIR, "fig1_rmse_vs_runtime.pdf"), fig)
end

function fig2_budget()
    # The final Figure 2 PDF is a retained frozen graphical artifact. Its
    # numerical source is the internal benchmark summary in `CHECKPOINT`.
    src = joinpath(@__DIR__, "inputs", "fig2_rmse_vs_budget_v3.pdf")
    cp(src, joinpath(FIGDIR, "fig2_rmse_vs_budget.pdf"); force=true)
end

function fig3_rules()
    df = CSV.read(INTEGRATION_SUMMARY, DataFrame)
    df = df[df.complete_paired_comparison .== true, :]
    cases = unique(df.case_id)
    rules = ["iid_mc", "halton", "sobol_owen", "lattice"]
    fig = Figure(size = (720, 430))
    ax = Axis(fig[1,1], xscale = log10, yscale = log10,
        xlabel = "median total runtime (s)", ylabel = "RMSE",
        title = "Integration-rule sensitivity")
    for rule in rules
        sub = df[df.rule .== rule, :]
        isempty(sub) && continue
        red = sub[sub.method .== "reduced_conditional", :]
        cur = sub[sub.method .== "current_indicator", :]
        scatter!(ax, Float64.(cur.median_total_runtime_s), Float64.(cur.rmse);
            color = (:gray35, 0.55), marker = :circle, markersize = 6,
            label = rule == first(rules) ? "Indicator cells" : nothing)
        scatter!(ax, Float64.(red.median_total_runtime_s), Float64.(red.rmse);
            color = (:dodgerblue3, 0.70), marker = :diamond, markersize = 7,
            label = rule == first(rules) ? "Reduced cells" : nothing)
    end
    axislegend(ax, position = :rt)
    Label(fig[2,1], "Rules: iid MC, randomized Halton, scrambled Sobol/Owen, shifted lattice.",
        tellwidth = false, fontsize = 10)
    save(joinpath(FIGDIR, "fig3_integration_rules.pdf"), fig)
end

function fig4_highdim()
    df = CSV.read(HIGHDIM_SUMMARY, DataFrame)
    by = Dict{Tuple{String,String,Int,String}, Dict{String,DataFrameRow}}()
    for r in eachrow(df)
        key = (String(r.case_id), String(r.event), Int(r.N), String(r.rule))
        get!(by, key, Dict{String,DataFrameRow}())[String(r.method)] = r
    end
    rows = NamedTuple[]
    for (_, m) in by
        haskey(m, "current_indicator") && haskey(m, "reduced_conditional") || continue
        c, red = m["current_indicator"], m["reduced_conditional"]
        push!(rows, (
            d = Int(red.d),
            rho = Float64(red.rho),
            event = String(red.event),
            N = Int(red.N),
            rmse_ratio = Float64(red.rmse) / Float64(c.rmse),
            runtime_ratio = Float64(red.median_total_runtime_s) / Float64(c.median_total_runtime_s),
        ))
    end
    out = DataFrame(rows)
    fig = Figure(size = (720, 360))
    ax1 = Axis(fig[1,1], yscale = log10, xlabel = "dimension d",
        ylabel = "Reduced / Indicator RMSE", title = "Accuracy ratio")
    ax2 = Axis(fig[1,2], yscale = log10, xlabel = "dimension d",
        ylabel = "Reduced / Indicator runtime", title = "Runtime ratio")
    for event in unique(out.event), rho in unique(out.rho)
        sub = sort(out[(out.event .== event) .&& (out.rho .== rho) .&& (out.N .== 256), :], :d)
        isempty(sub) && continue
        label = "$(event), ρ=$(rho)"
        marker = event == "cdf" ? :circle : :rect
        color = rho < 0.5 ? :dodgerblue3 : :firebrick3
        lines!(ax1, sub.d, sub.rmse_ratio; color, linewidth = 2, label)
        scatter!(ax1, sub.d, sub.rmse_ratio; color, marker, markersize = 8)
        lines!(ax2, sub.d, sub.runtime_ratio; color, linewidth = 2, label)
        scatter!(ax2, sub.d, sub.runtime_ratio; color, marker, markersize = 8)
    end
    hlines!(ax1, [1.0], color = :black, linestyle = :dash)
    hlines!(ax2, [1.0], color = :black, linestyle = :dash)
    axislegend(ax1, position = :rb, unique = true, rowgap = 2)
    save(joinpath(FIGDIR, "fig4_highdim_scaling.pdf"), fig)
end

function fig5_financial()
    df = CSV.read(FINANCIAL_SUMMARY, DataFrame)
    fig = Figure(size = (720, 360))
    events = ["joint_lower_10pct", "central_25_75_rectangle"]
    titles = Dict(
        "joint_lower_10pct" => "Lower 10% joint CDF",
        "central_25_75_rectangle" => "Central 25--75% rectangle",
    )
    for (j, event) in enumerate(events)
        sub = df[df.event .== event, :]
        ax = Axis(fig[1,j], xscale = log10, yscale = log10,
            xlabel = "median runtime (s)", ylabel = "RQMC standard error",
            title = titles[event])
        for method in ["current_indicator", "reduced_conditional"]
            m = sort(sub[sub.method .== method, :], :median_runtime_s)
            lines!(ax, Float64.(m.median_runtime_s), Float64.(m.rqmc_se);
                color = METHOD_COLOR[method], linewidth = 2, label = METHOD_LABEL[method])
            scatter!(ax, Float64.(m.median_runtime_s), Float64.(m.rqmc_se);
                color = METHOD_COLOR[method], marker = METHOD_MARKER[method], markersize = 8)
        end
    end
    axislegend(fig.content[1], position = :lb, unique = true)
    save(joinpath(FIGDIR, "fig5_financial_application.pdf"), fig)
end

fig1_internal_rmse_runtime()
fig2_budget()
fig3_rules()
fig4_highdim()
fig5_financial()

println("Wrote final manuscript figures to ", FIGDIR)
