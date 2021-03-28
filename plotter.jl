using GLMakie

function plot_agg(agglog::String; suffixes = ["energy", "sens", "repr", "speed"])
    df = CSV.File(agglog) |> DataFrame
    nspecies = maximum(df.species)
    gdf = groupby(df, :species)
    fig = Figure(resolution = (1600, 1000))
    cols = cgrad(:lighttest; alpha=0.4)
    

    plotline(ax, steps, data, color) = lines!(ax, steps, data; color, linewidth = 2)

    ax = fig[1, 1] = Axis(fig; ylabel = "Bact population")
    lines = [plotline(ax, gdf[i].step, gdf[i].nbact, cols[i/nspecies]) for i in 1:nspecies]

    for suffix in suffixes
        axm = fig[end+1, 1] = Axis(fig; ylabel = "Mean $suffix")
        axs = fig[end, 2] = Axis(fig; ylabel = "Std $suffix")

        for i in 1:nspecies
            plotline(axm, gdf[i].step, gdf[i][!, "μ_$suffix"], cols[i/nspecies])
            plotline(axs, gdf[i].step, gdf[i][!, "σ_$suffix"], cols[i/nspecies])
        end
    end
    fig
end
