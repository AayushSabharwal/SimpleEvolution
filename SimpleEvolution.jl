using Agents
using Random

const AXIAL_DIRECTIONS = [-1, 0, 1]

@agent Bacterium GridAgent{2} begin
    age::Int        # age, in iterations. They only live so long
    energy::Float64 # energy. Less than zero, it dies. More than threshold, it doubles
    sensory_radius::Float64 # how far it can see food
    reproduction_threshold::Float64 # see above
    speed::Int      # how far it can move in one iteration
    food_target::Dims{2}    # the food it can see
end

@agent Food GridAgent{2} begin
    food_cap::Float64   # maximum food value this has
    current_food::Float64   # current food value this has
    regen_rate::Float64 # how fast it regenerates. Faster regeneration means less resources going into multiplying
end

function initialize_model(
    dims::Dims{2},  # size of the space
    food_distribution::Array{Float64,2};    # food at each cell
    n_bacteria::Int = 100,  # number of bacteria
    energy::Float64 = 20.0, # initial energy of each bacteria
    sensory_radius::Float64 = 4.0,  # initial sensory radius of each bacteria
    threshold::Float64 = 50.0,  # initial reproduction threshold of each bacteria
    speed::Int = 1, # initial speed of each bacteria
    food_cap::Float64 = maximum(food_distribution), # initial food cap for all food
    regen_rate::Float64 = 10.0, # initial regen rate for all food
    lifetime::Int = 200,    # initial lifetime for all bacteria
    eat_rate_factor::Float64 = 5.0, # how fast all bacteria eat
    sensory_radius_cost_factor::Float64 = 1.0,  # energy cost multiplier for sensory radius
    distance_cost_factor::Float64 = 1.0,    # energy cost multipler for moving unit distance
    sensory_radius_std::Float64 = 1.0,  # random variation in inheriting sensory radius
    reproduction_threshold_std::Float64 = 1.0,  # random variation in inheriting reproduction threshold
    speed_std::Float64 = 2.0,   # random variation in inheriting speed
    reproduction_energy_cost::Float64 = 5.0,    # energy cost of reproduction
)
    # sanity checks
    @assert size(food_distribution) == dims
    @assert food_cap >= maximum(food_distribution)
    @assert lifetime > 0

    # reproducibility matters
    rng = MersenneTwister(42)

    # food lives in a separate model, since this improves performance
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
        :speed_std => speed_std,
    )

    # non-periodic euclidean space
    space = GridSpace(dims; periodic = false, metric = :euclidean)
    model = ABM(Bacterium, space; properties, rng)

    # add bacteria at random places
    for i ∈ 1:n_bacteria
        add_agent!(Bacterium, model, 0, energy, sensory_radius, threshold, speed, (-1, -1))
    end

    # add food wherever it is nonzero
    for pos ∈ eachindex(food_distribution)
        food_distribution[pos] > 0 || continue

        add_agent!(Food, model.food, food_cap, food_distribution[pos], regen_rate)
    end

    return model
end

"""
    eat(bact::Bacterium, food::Food, erf::Float64)

`erf` is `eat_rate_factor`. Makes `bact` eat from `food`.
"""
function eat(bact::Bacterium, food::Food, erf::Float64)
    eaten = min(erf * bact.speed, food.current_food)
    food.current_food -= eaten
    bact.energy += eaten
end

"""
    inherit(parameter::Float64, rng::MersenneTwister, std::Float64)
    inherit(parameter::Int, rng::MersenneTwister, std::Float64)

Simulates inheritance of `parameter` with given `standard deviation` using normal distribution.
"""
inherit(parameter::Float64, rng::MersenneTwister, std::Float64) =
    max(parameter + randn(rng) * std, 1.0)
inherit(parameter::Int, rng::MersenneTwister, std::Float64) =
    max(parameter + floor(Int, randn(rng) * std), 1)

"""
    reproduce(bact::Bacterium, model::ABM)

Simulates reproduction of bacteria. Energy is halved after subtracting a cost. Inherited features are
varied.
"""
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
    # die, of age or startvation
    (bact.age > model.lifetime || bact.energy <= 0) && (kill_agent!(bact, model); return)

    # grow older
    bact.age += 1

    # reproduce if possible
    bact.energy >= bact.reproduction_threshold && reproduce(bact, model)

    # if currently on a food source
    if !isempty(bact.pos, model.food)
        # get the food
        food = first(agents_in_position(bact.pos, model.food))
        # if there's no food on it
        if food.current_food == 0
            # can't eat, so look for something else
            bact.food_target = (-1, -1)
        else
            # eat
            eat(bact, food, model.eat_rate_factor)
            return
        end
    end

    # if bact doesn't see any food
    if bact.food_target == (-1, -1)
        best_id = -1
        # look for the best food nearby
        for id in nearby_ids(bact.pos, model.food, bact.sensory_radius)
            (
                best_id == -1 ||
                model.food[best_id].current_food < model.food[id].current_food
            ) || continue
            best_id = id
        end
        best_id != -1 && (bact.food_target = (-1, -1))
    end

    # if we found some food or already had eyes on it
    if bact.food_target != (-1, -1)
        # move toward target
        delta = bact.food_target .- bact.pos
        #move as many steps diagonally as possible
        diag = min(abs.(delta)..., bact.speed)

        steps = diag

        movement = sign.(delta) .* diag
        # move the rest steps axially, if possible
        delta = delta .- movement
        while steps < bact.speed && any(delta .> 0)
            movement = movement .+ sign.(delta)
            delta = delta .- sign.(delta)
            steps += 1
        end
    else
        # move randomly if there's no food nearby
        movement = rand.(model.rng, (AXIAL_DIRECTIONS, AXIAL_DIRECTIONS)) .* bact.speed
    end

    # move
    move_agent!(bact, clamp.(bact.pos .+ movement, 1, size(model.space.s)), model)

    # subtract cost of movement and looking
    bact.energy -=
        model.sensory_radius_cost_factor * bact.sensory_radius +
        model.distance_cost_factor * max(abs.(movement)...)
end
