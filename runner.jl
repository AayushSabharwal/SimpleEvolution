using DataFrames
using StatsBase
using CSV
using ProgressMeter

function run!(
    model::ABM,
    agglog::String,
    bactlog::String,
    foodlog::String;
    nsteps::Int = 2000,
    log_period::Int = 1,
    chunk_size::Int = 1000,
)
    agg_data = DataFrame(
        step = Int[],
        nbact = Int[],
        nfood = Int[],
        tot_food = Float64[],
        μ_energy = Float64[],
        σ_energy = Float64[],
        μ_sens = Float64[],
        σ_sens = Float64[],
        μ_repr = Float64[],
        σ_repr = Float64[],
        μ_speed = Float64[],
        σ_speed = Float64[],
        μ_food_cap = Float64[],
        σ_food_cap = Float64[],
        μ_regen_rate = Float64[],
        σ_regen_rate = Float64[],
    )
    bactlogdf = DataFrame(
        step = Int[],
        age = Int[],
        sens = Float64[],
        repr = Float64[],
        speed = Float64[],
    )
    foodlogdf = DataFrame(
        step = Int[],
        food_cap = Float64[],
        regen_rate = Float64[],
        x = Int[],
        y = Int[],
    )
    appendbact = false
    appendfood = false
    appendagg = false

    p = Progress(nsteps; barglyphs = BarGlyphs("[=> ]"), color = :blue)
    for i = 0:nsteps
        push!(
            agg_data,
            (
                i,
                nagents(model),
                nagents(model.food),
                sum(a.current_food for a in allagents(model.food)),
                mean_and_std([a.energy for a in allagents(model)])...,
                mean_and_std([a.sensory_radius for a in allagents(model)])...,
                mean_and_std([a.reproduction_threshold for a in allagents(model)])...,
                mean_and_std([a.speed for a in allagents(model)])...,
                mean_and_std([a.food_cap for a in allagents(model.food)])...,
                mean_and_std([a.regen_rate for a in allagents(model.food)])...,
            ),
        )
        if size(agg_data, 1) >= chunk_size
            CSV.write(agglog, agg_data; append = appendagg)
            appendagg || (appendagg = true)
            empty!(agg_data)
        end

        if i % log_period == 0
            for a in allagents(model)
                push!(
                    bactlogdf,
                    (i, a.age, a.sensory_radius, a.reproduction_threshold, a.speed),
                )
            end
            for a in allagents(model.food)
                push!(foodlogdf, (i, a.food_cap, a.regen_rate, a.pos[1], a.pos[2]))
            end
            if size(bactlogdf, 1) >= chunk_size
                CSV.write(bactlog, bactlogdf; append = appendbact)
                appendbact || (appendbact = true)
                empty!(bactlogdf)
            end
            if size(foodlogdf, 1) >= chunk_size
                CSV.write(foodlog, foodlogdf; append = appendfood)
                appendfood || (appendfood = true)
                empty!(foodlogdf)
            end
        end
        Agents.step!(model, agent_step!, food_step!)
        ProgressMeter.next!(p)
    end
    size(agg_data, 1) > 0 && CSV.write(agglog, agg_data; append = appendagg)
    size(bactlogdf, 1) > 0 && CSV.write(bactlog, bactlogdf; append = appendbact)
    size(foodlogdf, 1) > 0 && CSV.write(foodlog, foodlogdf; append = appendfood)
    nothing
end
