using CairoMakie
using GLMakie
using Colors

# Theme used
const BG = colorant"#191919"
const AX = colorant"#ababac"
const TX = colorant"#cdcdcf"
const GR = colorant"#717280"
const THEME = Theme(
    backgroundcolor = BG,
    Axis = (
        backgroundcolor = :transparent,
        bottomspinecolor = AX,
        leftspinecolor = AX,
        rightspinecolor = AX,
        topspinecolor = AX,
        xgridcolor = GR,
        xtickcolor = AX,
        xlabelcolor = TX,
        xticklabelcolor = TX,
        ygridcolor = GR,
        ytickcolor = AX,
        ylabelcolor = TX,
        yticklabelcolor = TX,
    ),
    Legend = (bgcolor = :transparent, framecolor = AX, labelcolor = TX),
    Label = (color = TX,),
)
set_theme!(THEME)
# Full names of parameters
const PARAM_MAP =
    Dict("repr" => "Reproduction Threshold", "speed" => "Speed", "sens" => "Sensory Radius")

# Plot the food as a heatmap, with a slider for interactivity
function plot_food(foodlog::String)
    fig = Figure(resolution = (600, 600))
    food = h5open(foodlog, "r")
    maxtime = size(food["log"], 1)
    time = Slider(fig[1, 1], range = 1:maxtime, startvalue = 1)
    cfood = @lift(food["log"][$(time.value), :, :])
    ax = fig[0, 1] = Axis(fig)
    heatmap!(ax, cfood, colorrange = (0, food["cap"][]))

    fig, food
end

# Record the food map playing through steps
function record_food(
    foodlog::String,
    filename::String;
    duration::Real = 10,
    nsteps::Int = -1,
)
    CairoMakie.activate!()
    fig = Figure(resolution = (600, 600), figure_padding = 0)
    food = h5open(foodlog, "r")
    maxtime = size(food["log"], 1)
    nsteps = nsteps <= 0 ? maxtime : min(nsteps, maxtime)
    ax = fig[1, 1] = Axis(fig; aspect = AxisAspect(1))
    hidedecorations!(ax)
    time = Observable(1)
    cfood = @lift(food["log"][$(time), :, :])
    titletext = @lift("Step: $($(time))")
    fig[0, 1] = Label(fig, titletext; tellwidth = false)
    heatmap!(ax, cfood, colorrange = (0, food["cap"][]))
    record(fig, filename, 1:nsteps; framerate = floor(Int, nsteps / duration)) do i
        time[] = i
    end
    close(food)
    GLMakie.activate!()
end

# Plot data of the simulation
function plot_agent_scatter(
    dir::String;
    params = ["sens", "repr", "speed"],
    fig = Figure(resolution = (1600, 900)),
    plotbact = false,
)
    agg_df = CSV.File(joinpath(dir, "agg.csv")) |> DataFrame
    n_species = maximum(agg_df.species)
    stepcount = maximum(agg_df.step)
    print(agg_df.step)
    print(stepcount)
    print(typeof(stepcount))
    agg_df = groupby(agg_df, :species; sort = true)

    plotline(ax, steps, data, color) = lines!(ax, steps, data; color, linewidth = 2)

    if plotbact
        df =
            CSV.File(joinpath(dir, "bact.csv"); select = ["step", "species", params...]) |>
            DataFrame
        df = groupby(df, :species; sort = true)
        scatter_cols = distinguishable_colors(n_species, HSV(0, 0.6, 1))
    end

    sp_cols = distinguishable_colors(n_species, HSV(0, 0.6, 1))

    ax =
        fig[1, :] = Axis(
            fig;
            ylabel = "Population",
            xticks = 0:500:stepcount,
            xticksvisible = false,
            xticklabelsvisible = false,
        )
    xlims!(ax, (0, stepcount))
    lines = [plotline(ax, agg_df[i].step, agg_df[i].nbact, sp_cols[i]) for i = 1:n_species]

    for (ind, par) in enumerate(params)
        ax =
            fig[ind+1, :] = Axis(
                fig;
                ylabel = PARAM_MAP[par],
                xticks = 0:500:stepcount,
                xticksvisible = ind == length(params),
                xticklabelsvisible = ind == length(params),
            )
        for i = 1:n_species
            xlims!(ax, (0, stepcount))
            plotbact && scatter!(
                ax,
                df[i].step,
                df[i][!, par];
                color = (scatter_cols[i], 0.01),
                markersize = 4,
                strokewidth = 0,
            )
            plotline(ax, agg_df[i].step, agg_df[i][!, "μ_$par"], sp_cols[i])
            band!(
                ax,
                agg_df[i].step,
                agg_df[i][!, "μ_$par"] .- agg_df[i][!, "σ_$par"],
                agg_df[i][!, "μ_$par"] .+ agg_df[i][!, "σ_$par"],
                color = (sp_cols[i], 0.25),
            )
        end
    end

    leg =
        fig[length(params)+2, :] = Legend(
            fig,
            lines,
            ["Mean" for _ = 1:n_species];
            orientation = :horizontal,
            tellwidth = false,
            tellheight = true,
        )
    fig
end

# Multiple `plot_agent_scatter` side by side
function comparative_agent_scatter(dirs::Vector{String}; kwargs...)
    fig = Figure(resolution = (1200, 1200))

    for (i, dir) in enumerate(dirs)
        plot_agent_scatter(dir; fig = fig[1, i], kwargs...)
    end
    fig
end
