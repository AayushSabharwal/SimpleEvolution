using Agents

@agent Bacterium GridAgent{2} begin
    age::Int
    energy::Float64
    sensory_radius::Float64
    reproduction_threshold::Float64
    speed::Int
    food_target::Dims{2}
end

@agent Food GridAgent{2} begin
    food_cap::Float64
    current_food::Float64
    regen_rate::Float64
end

function initialize_model(
    dims::Dims{2},
    food_distribution::Array{Float64,2};
    n_bacteria::Int = 100,
    energy::Float64 = 20.0,
    sensory_radius::Float64 = 4.0,
    threshold::Float64 = 50.0,
    speed::Int = 1,
    food_cap::Float64 = maximum(food_distribution),
    regen_rate::Float64 = 10.0,
    lifetime::Int = 200,
    eat_rate_factor::Float64 = 5.0,
    sensory_radius_cost_factor::Float64 = 1.0,
    distance_cost_factor::Float64 = 1.0,
    sensory_radius_std::Float64 = 1.0,
    reproduction_threshold_std::Float64 = 1.0,
    speed_std::Float64 = 2.0,
    reproduction_energy_cost::Float64 = 5.0,
)
    @assert size(food_distribution) == dims
    @assert food_cap >= maximum(food_distribution)
    @assert lifetime > 0

    rng = MersenneTwister(42)

    food_space = GridSpace(dims; periodic = false)
    food_model = ABM(Food, food_space; rng)

    properties = Dict(
        :eat_rate_factor => eat_rate_factor,
        :lifetime => lifetime,
        :food => food_model,
        :sensory_radius_cost_factor => sensory_radius_cost_factor,
        :distance_cost_factor => distance_cost_factor,
        :sensory_radius_std => sensory_radius_std,
        :reproduction_threshold_std => reproduction_threshold_std,
        :reproduction_energy_cost => reproduction_energy_cost,
    )

    space = GridSpace(dims; periodic = false, metric = :euclidean)
    model = ABM(Bacterium, space; properties, rng)


    for i ∈ 1:n_bacteria
        add_agent!(Bacterium, model, 0, energy, sensory_radius, threshold, speed, (-1, -1))
    end

    for pos ∈ eachindex(food_distribution)
        food_distribution[pos] > 0 || continue

        add_agent!(Food, model.food, food_cap, food_distribution[pos], regen_rate)
    end

    return model
end

function eat(bact::Bacterium, food::Food, erf::Float64)
    eaten = min(erf * bact.speed, food.current_food)
    food.current_food -= eaten
    bact.energy += eaten
end

inherit(parameter::Float64, rng::MersenneTwister, std::Float64) =
    max(parameter + randn(rng) * std, 1.0)
inherit(parameter::Int, rng::MersenneTwister, std::Float64) =
    max(parameter + floor(Int, randn(rng) * std), 1)

function reproduce(bact::Bacterium, model::ABM)
    add_agent!(
        bact.pos,
        Bacterium,
        model,
        0,
        (bact.energy - model.reproduction_energy_cost) / 2.0,
        inherit(bact.sensory_radius, model.rng, model.sensory_radius_std),
        inherit(bact.reproduction_threshold, model.rng, model.reproduction_threshold_std),
        inherit(bact.speed, model.rng, model.speed_std),
        bact.food_target,
    )

    bact.age = 0
    bact.energy -= model.reproduction_energy_cost
    bact.energy /= 2.0

    bact.sensory_radius = inherit(bact.sensory_radius, model.rng, model.sensory_radius_std)
    bact.reproduction_threshold =
        inherit(bact.reproduction_threshold, model.rng, model.reproduction_threshold_std)
    bact.speed = inherit(bact.speed_std, model.rng, model.speed_std)
end

function agent_step!(bact::Bacterium, model::ABM)
    bact.age < model.lifetime || (kill_agent!(bact, model); return)

    bact.energy >= bact.reproduction_threshold && reproduce(bact, model)

    isempty(bact.pos, model.food) || (
        eat(bact, first(agents_in_position(bact.pos, model.food)), model.eat_rate_factor); return
    )

    bact.age += 1
    if bact.food_target == (-1, -1)
        best_id = -1
        for id in nearby_ids(bact.pos, model.food, bact.sensory_radius)
            (
                best_id == -1 ||
                model.food[best_id].current_food < model.food[id].current_food
            ) || continue
            best_id = id
        end
        best_id != -1 && (bact.food_target = (-1, -1))
    end

    if bact.food_target != (-1, -1)
        delta = bact.food_target .- bact.pos
        diag = min(abs.(delta)..., bact.speed)

        steps = diag

        movement = sign.(delta) .* diag
        delta = delta .- movement
        while steps < bact.speed && any(delta .> 0)
            movement = movement .+ sign.(delta)
            delta = delta .- sign.(delta)
            steps += 1
        end

        move_agent!(bact, model, bact.pos .+ movement)
    else
        movement = rand.(model.rng, ([-1, 0, 1], [-1, 0, 1])) .* bact.speed
    end

    bact.energy -=
        model.sensory_radius_cost_factor * bact.sensory_radius +
        model.distance_cost_factor * max(abs.(movement)...)
    bact.energy <= 0 && kill_agent!(bact, model)
end
