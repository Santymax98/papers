using Statistics

const OUTDIR = @__DIR__
const SUMMARY = joinpath(OUTDIR, "integration_rule_summary.csv")
const FIGDIR = joinpath(OUTDIR, "figures")
mkpath(FIGDIR)

function csv_parse_line(line::AbstractString)
    fields = String[]
    buf = IOBuffer()
    inq = false
    i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if inq
            if c == '"'
                ni = nextind(line, i)
                if ni <= lastindex(line) && line[ni] == '"'
                    write(buf, '"')
                    i = ni
                else
                    inq = false
                end
            else
                write(buf, c)
            end
        elseif c == '"'
            inq = true
        elseif c == ','
            push!(fields, String(take!(buf)))
        else
            write(buf, c)
        end
        i = nextind(line, i)
    end
    push!(fields, String(take!(buf)))
    return fields
end

function read_csv(path)
    lines = readlines(path)
    header = csv_parse_line(first(lines))
    rows = Dict{String,String}[]
    for line in Iterators.drop(lines, 1)
        isempty(strip(line)) && continue
        push!(rows, Dict(zip(header, csv_parse_line(line))))
    end
    return rows
end

function esc(s)
    replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")
end

function write_scatter_svg(path, rows; xfield, yfield, title, xlabel, ylabel)
    rows = filter(r -> parse(Bool, r["complete_paired_comparison"]) &&
        r["method"] == "reduced_conditional", rows)
    currows = filter(r -> parse(Bool, r["complete_paired_comparison"]) &&
        r["method"] == "current_indicator", read_csv(SUMMARY))
    allpts = Tuple{Float64,Float64,String,String,String}[]
    colors = Dict("iid_mc"=>"#666666", "randomized_halton"=>"#0072B2",
                  "sobol_owen"=>"#009E73", "shifted_lattice"=>"#D55E00")
    symbols = Dict("current_indicator"=>"open", "reduced_conditional"=>"filled")
    for r in vcat(currows, rows)
        x = log10(parse(Float64, r[xfield]))
        y = log10(parse(Float64, r[yfield]))
        push!(allpts, (x, y, r["rule"], r["method"], r["case_id"]))
    end
    xs = first.(allpts)
    ys = getindex.(allpts, 2)
    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)
    padx = max(0.05, 0.06 * (xmax - xmin))
    pady = max(0.05, 0.06 * (ymax - ymin))
    xmin -= padx; xmax += padx; ymin -= pady; ymax += pady
    W, H = 920, 620
    l, r, t, b = 82, 24, 48, 82
    sx(x) = l + (x - xmin) / (xmax - xmin) * (W - l - r)
    sy(y) = H - b - (y - ymin) / (ymax - ymin) * (H - t - b)
    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$W" height="$H" viewBox="0 0 $W $H">""")
        println(io, """<rect width="100%" height="100%" fill="white"/>""")
        println(io, """<text x="$(W/2)" y="28" text-anchor="middle" font-family="Helvetica" font-size="18">$(esc(title))</text>""")
        println(io, """<line x1="$l" y1="$(H-b)" x2="$(W-r)" y2="$(H-b)" stroke="black"/>""")
        println(io, """<line x1="$l" y1="$t" x2="$l" y2="$(H-b)" stroke="black"/>""")
        println(io, """<text x="$(W/2)" y="$(H-25)" text-anchor="middle" font-family="Helvetica" font-size="14">$(esc(xlabel))</text>""")
        println(io, """<text x="22" y="$(H/2)" transform="rotate(-90 22 $(H/2))" text-anchor="middle" font-family="Helvetica" font-size="14">$(esc(ylabel))</text>""")
        for (x, y, rule, method, case_id) in allpts
            fill = method == "current_indicator" ? "white" : colors[rule]
            stroke = colors[rule]
            radius = method == "current_indicator" ? 4 : 3
            println(io, """<circle cx="$(sx(x))" cy="$(sy(y))" r="$radius" fill="$fill" stroke="$stroke" stroke-width="1.4"><title>$(esc(case_id)) $(esc(rule)) $(esc(method))</title></circle>""")
        end
        lx, ly = W - 250, 70
        for (i, rule) in enumerate(sort(collect(keys(colors))))
            y = ly + 22i
            println(io, """<circle cx="$lx" cy="$y" r="4" fill="$(colors[rule])" stroke="$(colors[rule])"/>""")
            println(io, """<text x="$(lx+14)" y="$(y+5)" font-family="Helvetica" font-size="12">$(esc(rule))</text>""")
        end
        println(io, """<circle cx="$lx" cy="$(ly+130)" r="4" fill="white" stroke="black"/>""")
        println(io, """<text x="$(lx+14)" y="$(ly+135)" font-family="Helvetica" font-size="12">current indicator</text>""")
        println(io, """<circle cx="$lx" cy="$(ly+152)" r="4" fill="black" stroke="black"/>""")
        println(io, """<text x="$(lx+14)" y="$(ly+157)" font-family="Helvetica" font-size="12">reduced conditional</text>""")
        println(io, "</svg>")
    end
end

rows = read_csv(SUMMARY)
write_scatter_svg(
    joinpath(FIGDIR, "integration_rule_rmse_vs_budget.svg"),
    rows;
    xfield="N",
    yfield="rmse",
    title="Phase 3: RMSE versus budget",
    xlabel="log10(N)",
    ylabel="log10(RMSE)",
)
write_scatter_svg(
    joinpath(FIGDIR, "integration_rule_rmse_vs_runtime.svg"),
    rows;
    xfield="median_total_runtime_s",
    yfield="rmse",
    title="Phase 3: RMSE versus total runtime",
    xlabel="log10(total runtime seconds)",
    ylabel="log10(RMSE)",
)
println("Wrote Phase 3 SVG figures to $FIGDIR")

