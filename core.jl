const AXIAL_DIRECTIONS = [-1, 0, 1]
const NEIGHBORHOOD = [a for a in Iterators.product([-1:1 for φ = 1:2]...) if a != (0, 0)]

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

"""
    inherit(parameter::Float64, rng::MersenneTwister, std::Float64)
    inherit(parameter::Int, rng::MersenneTwister, std::Float64)

Simulates inheritance of `parameter` with given `standard deviation` using normal distribution.
"""
inherit(parameter::Float64, rng::MersenneTwister, std::Float64) =
    max(parameter + randn(rng) * std, 1.0)
inherit(parameter::Int, rng::MersenneTwister, std::Float64) =
    max(parameter + floor(Int, randn(rng) * std), 1)

function initialize_model(
    dims::Dims{2},  # size of the space
    food_distribution::Array{Float64,2};    # food at each cell
    n_bacteria::Int = 100,  # number of bacteria
    initial_bacterium::Bacterium = Bacterium(-1, (-1, -1), 0, 20.0, 4.0, 50.0, 1, (-1, -1)),
    initial_food::Food = Food(-1, (-1, -1), maximum(food_distribution), 0.0, 10.0),
    lifetime::Int = 200,    # initial lifetime for all bacteria
    eat_rate_factor::Float64 = 5.0, # how fast all bacteria eat
    sensory_radius_cost_factor::Float64 = 1.0,  # energy cost multiplier for sensory radius
    distance_cost_factor::Float64 = 1.0,    # energy cost multipler for moving unit distance
    sensory_radius_std::Float64 = 1.0,  # random variation in inheriting sensory radius
    reproduction_threshold_std::Float64 = 1.0,  # random variation in inheriting reproduction threshold
    speed_std::Float64 = 2.0,   # random variation in inheriting speed
    reproduction_energy_cost::Float64 = 5.0,    # energy cost of reproduction
    food_spread_coefficient::Float64 = 3.0,
    food_reproduction_cost::Float64 = 10.0,
    food_cap_std::Float64 = 5.0,
    regen_rate_std::Float64 = 5.0,
)
    # sanity checks
    @assert size(food_distribution) == dims
    @assert initial_food.food_cap >= maximum(food_distribution)
    @assert lifetime > 0

    # reproducibility matters
    rng = MersenneTwister(42)

    food_properties = Dict(
        :spread_coefficient => food_spread_coefficient,
        :reproduction_cost => food_reproduction_cost,
        :food_cap_std => food_cap_std,
        :regen_rate_std => regen_rate_std,
    )

    # food lives in a separate model, since this improves performance
    food_space = GridSpace(dims; periodic = false)
    food_model = ABM(Food, food_space; properties = food_properties, rng)

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
        add_agent!(
            Bacterium,
            model,
            initial_bacterium.age,
            initial_bacterium.energy,
            initial_bacterium.sensory_radius,
            initial_bacterium.reproduction_threshold,
            initial_bacterium.speed,
            (-1, -1),
        )
    end

    # add food wherever it is nonzero
    for pos ∈ eachindex(food_distribution)
        food_distribution[pos] > 0 || continue

        add_agent!(
            Food,
            model.food,
            initial_food.food_cap,
            food_distribution[pos],
            initial_food.regen_rate,
        )
    end

    return model
end
