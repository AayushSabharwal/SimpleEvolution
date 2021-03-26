using GLMakie

function plot_agg(agglog::String)
    df = CSV.File(agglog) |> DataFrame
    fig = Figure(resolution = (1600, 1000))

    plotline(ax, data, color) = lines!(ax, df.step, data; color, linewidth = 2)

    ax = fig[1, 1] = Axis(fig; ylabel = "Bact population")
    redl = plotline(ax, df.nbact, (:red, 0.4))

    ax = fig[1, 2] = Axis(fig; ylabel = "Total food")
    redl = plotline(ax, df.tot_food, (:blue, 0.4))

    ax = fig[end+1, 1] = Axis(fig; ylabel = "Energy")
    plotline(ax, df.μ_energy, (:red, 0.4))
    plotline(ax, df.σ_energy, (:blue, 0.4))

    ax = fig[end, 2] = Axis(fig; ylabel = "Sensory Radius")
    plotline(ax, df.μ_sens, (:red, 0.4))
    plotline(ax, df.σ_sens, (:blue, 0.4))

    ax = fig[end+1, 1] = Axis(fig; ylabel = "Repr. Threshold")
    plotline(ax, df.μ_repr, (:red, 0.4))
    plotline(ax, df.σ_repr, (:blue, 0.4))

    ax = fig[end, 2] = Axis(fig, ylabel = "Speed")
    plotline(ax, df.μ_speed, (:red, 0.4))
    plotline(ax, df.σ_speed, (:blue, 0.4))

    ax = fig[end+1, 1] = Axis(fig, ylabel = "Food Cap")
    plotline(ax, df.μ_food_cap, (:red, 0.4))
    plotline(ax, df.σ_food_cap, (:blue, 0.4))

    ax = fig[end, 2] = Axis(fig, ylabel = "Regen Rate")
    redl = plotline(ax, df.μ_regen_rate, (:red, 0.4))
    bluel = plotline(ax, df.σ_regen_rate, (:blue, 0.4))

    fig[end+1, 1:2] = Legend(
        fig,
        [redl, bluel],
        ["mean", "std"];
        orientation = :horizontal,
        tellheight = true,
    )
    fig
end
