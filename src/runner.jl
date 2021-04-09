using DataFrames
using StatsBase
using CSV
using ProgressMeter
using OnlineStats
using HDF5

function run!(
    model::ABM,
    agglog::String,
    bactlog::String,
    foodlog::String;
    nsteps::Int = 2000,
    log_period::Int = 1,
    log_bact::Bool = true,
    log_food::Bool = true,
    chunk_size::Int = 1000,
    float_precision::Int = 2,
)
    agg_data = DataFrame(
        step = Int[],
        species = Int[],
        nbact = Int[],
        μ_energy = Float64[],
        σ_energy = Float64[],
        μ_sens = Float64[],
        σ_sens = Float64[],
        μ_repr = Float64[],
        σ_repr = Float64[],
        μ_speed = Float64[],
        σ_speed = Float64[],
    )
    bactlogdf = DataFrame(
        step = Int[],
        age = Int[],
        species = Int[],
        energy = Float64[],
        sens = Float64[],
        repr = Float64[],
        speed = Float64[],
    )
    
    agglog = open(datadir(agglog), "w+")
    log_bact && (bactlog = open(datadir(bactlog), "w+"))
    foodlog = datadir(foodlog)

    if log_food
        isfile(foodlog) && rm(foodlog)
        foodf = h5open(foodlog, "w")
        foodf["log"] = zeros((floor(Int, nsteps/log_period)+1, size(model.food)...))
        fdata = HDF5.readmmap(foodf["log"])
    end
    appendbact = false
    appendagg = false

    p = Progress(nsteps; barglyphs = BarGlyphs("[=> ]"), color = :blue)
    for i = 0:nsteps
        count = [0 for _ in 1:model.nspecies]
        energy_means = [Mean() for _ in 1:model.nspecies]
        energy_vars = [Variance() for _ in 1:model.nspecies]
        sens_means = [Mean() for _ in 1:model.nspecies]
        sens_vars = [Variance() for _ in 1:model.nspecies]
        repr_means = [Mean() for _ in 1:model.nspecies]
        repr_vars = [Variance() for _ in 1:model.nspecies]
        speed_means = [Mean() for _ in 1:model.nspecies]
        speed_vars = [Variance() for _ in 1:model.nspecies]

        for a in allagents(model)
            count[a.species] += 1
            fit!(energy_means[a.species], a.energy)
            fit!(energy_vars[a.species], a.energy)
            fit!(sens_means[a.species], a.sensory_radius)
            fit!(sens_vars[a.species], a.sensory_radius)
            fit!(repr_means[a.species], a.reproduction_threshold)
            fit!(repr_vars[a.species], a.reproduction_threshold)
            fit!(speed_means[a.species], a.speed)
            fit!(speed_vars[a.species], a.speed)
        end
        for sp in 1:model.nspecies
            push!(
                agg_data,
                (
                    i,
                    sp,
                    count[sp],
                    round(value(energy_means[sp]); digits = float_precision),
                    round(√value(energy_vars[sp]); digits = float_precision),
                    round(value(sens_means[sp]); digits = float_precision),
                    round(√value(sens_vars[sp]); digits = float_precision),
                    round(value(repr_means[sp]); digits = float_precision),
                    round(√value(repr_vars[sp]); digits = float_precision),
                    round(value(speed_means[sp]); digits = float_precision),
                    round(√value(speed_vars[sp]); digits = float_precision),
                ),
            )
        end
        if size(agg_data, 1) >= chunk_size
            CSV.write(agglog, agg_data; append = appendagg)
            appendagg || (appendagg = true)
            empty!(agg_data)
        end

        if log_period > 0 && i % log_period == 0
            if log_bact
                for a in allagents(model)
                    push!(
                        bactlogdf,
                        (i, a.age, a.species, round(a.energy; digits = float_precision), round(a.sensory_radius; digits = float_precision), round(a.reproduction_threshold; digits = float_precision), round(a.speed; digits = float_precision),)
                    )
                end
                if size(bactlogdf, 1) >= chunk_size
                    CSV.write(bactlog, bactlogdf; append = appendbact)
                    appendbact || (appendbact = true)
                    empty!(bactlogdf)
                end
            end
            
            log_food && (fdata[floor(Int, i/log_period)+1,:,:] = model.food)
        end
        Agents.step!(model, agent_step!, food_step!)
        ProgressMeter.next!(p)
    end
    
    size(agg_data, 1) > 0 && CSV.write(agglog, agg_data; append = appendagg)
    size(bactlogdf, 1) > 0 && CSV.write(bactlog, bactlogdf; append = appendbact)

    if log_food
        foodf["cap"] = model.food_data.food_cap
        print("CCC")
        close(foodf)
    end
    nothing
end

# from: https://gist.github.com/oyd11/392bfdc4fabae3d7361b5e52c961edd8
f_helper(x) = x
f_helper(d::Dict) = Dict(Symbol(k) => f_helper(v) for (k, v) in d)
symbol_dict(d::Dict) = f_helper(d)

function runconfig(
    config::Dict,
    dirname::String,
    foodpath::String;
    nsteps::Int = 2000,
    log_period::Int = 10,
    log_bact::Bool = true,
    log_food::Bool = true,
    chunk_size::Int = 1000,
)
    mkpath(datadir(dirname))
    config = symbol_dict(config)
    println(config)
    if config[:initial_bacterium] isa Array
        if length(config[:initial_bacterium]) == 0
            config[:initial_bacterium] = Bacterium()
        else
            ibact = Bacterium[]
            for bconf in config[:initial_bacterium]
                bconf = symbol_dict(bconf)
                println(bconf...)
                push!(ibact, Bacterium(; bconf...))
            end
            config[:initial_bacterium] = ibact
        end
    elseif  config[:initial_bacterium] isa Dict
        config[inital_bacterium] = Bacterium(; config[:initial_bacterium]...)
    end
    
    config[:n_bacteria] isa Vector && (config[:n_bacteria] = Vector{Int}(config[:n_bacteria]))

    config[:food_data] = FoodData(; config[:food_data]...)
    
    food = image_to_foodmap(foodpath)

    model = initialize_model(food .* config[:food_data].food_cap; config...)

    run!(model, joinpath(dirname, "agg.csv"), joinpath(dirname, "bact.csv"), joinpath(dirname, "food.hdf5"); nsteps, log_period, log_bact, log_food, chunk_size)
end