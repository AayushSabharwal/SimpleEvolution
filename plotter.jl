using GLMakie

function plot_agg(agglog::String; suffixes = ["energy", "sens", "repr", "speed"])
    df = CSV.File(agglog) |> DataFrame
    nspecies = maximum(df.species)
    gdf = groupby(df, :species; sort = true)
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

function plot_food(foodlog::String)
    fig = Figure(resolution = (600, 600))
    food = h5open(foodlog, "r")
    maxtime = size(food["log"], 1)
    time = Slider(fig[1, 1], range = 1:maxtime, startvalue = 1)
    cfood = @lift(food["log"][$(time.value), :, :])
    ax = fig[0, 1] = Axis(fig)
    heatmap!(ax, cfood)
    
    fig, food
end

function plot_agent_scatter(bactlog::String, agglog::String; params = ["sens", "repr", "speed"])
    agg_df = CSV.File(agglog) |> DataFrame
    n_species = maximum(agg_df.species)
    agg_df = groupby(agg_df, :species; sort = true)
    cols = cgrad(:RdYlGn_4; alpha=0.8)

    plotline(ax, steps, data, color) = lines!(ax, steps, data; color, linewidth = 2)

    df = CSV.File(bactlog; select = ["step", "species", params...]) |> DataFrame
    df = groupby(df, :species; sort = true)
    
    sp_cols = cgrad(:Spectral_6; alpha=0.1)
    fig = Figure(resolution = (600, 600))
    
    l1 = l2 = nothing
    for (ind, par) in enumerate(params)
        for i in 1:n_species
            ax = fig[1+ind, i] = Axis(fig; ylabel = par)
            xlims!(ax, (0, maximum(df[1].step)))
            scatter!(ax, df[i].step, df[i][!, par]; color = sp_cols[i/n_species], markersize = 4, strokewidth = 0)
            l1 = plotline(ax, agg_df[i].step, agg_df[i][!, "μ_$par"], cols[0.])
            l2 = plotline(ax, agg_df[i].step, agg_df[i][!, "σ_$par"], cols[1.])
        end
    end

    ax = fig[1, :] = Axis(fig; ylabel = "Bact population")
    lines = [plotline(ax, agg_df[i].step, agg_df[i].nbact, (sp_cols[i/n_species], 1.0)) for i in 1:n_species]
    
    fig[end+1, :] = Legend(fig, [l1, l2], ["Mean", "Std"]; orientation = :horizontal)
    fig
end