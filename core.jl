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
    regen_rate::Float64 # how fast it regenerates
end

Bacterium(;
    energy::Float64 = 20.0,
    sensory_radius::Float64 = 4.0,
    reproduction_threshold::Float64 = 80.0,
    speed::Int = 3,
) = Bacterium(
    -1,
    (-1, -1),
    0,
    energy,
    sensory_radius,
    reproduction_threshold,
    speed,
    (-1, -1),
)

Bacterium(id::Int, pos::Dims{2}, template::Bacterium) = Bacterium(
    id,
    pos,
    0,
    template.energy,
    template.sensory_radius,
    template.reproduction_threshold,
    template.speed,
    (-1, -1),
)

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
    n_bacteria::Int = 10,  # number of bacteria
    initial_bacterium::Bacterium = Bacterium(),
    food_cap_multiplier::Float64 = 2.0,
    initial_regen_rate::Float64 = 80.0,
    # model properties
    lifetime::Int = 200,    # initial lifetime for all bacteria
    eat_rate_factor::Float64 = 5.0, # how fast all bacteria eat
    sensory_radius_cost::Float64 = 5.0,  # energy cost multiplier for sensory radius
    distance_cost::Float64 = 5.0,    # energy cost multipler for moving unit distance
    σ_sensory_radius::Float64 = 1.0,  # random variation in inheriting sensory radius
    σ_reproduction_threshold::Float64 = 1.0,  # random variation in inheriting reproduction threshold
    σ_speed::Float64 = 2.0,   # random variation in inheriting speed
    reproduction_energy_cost::Float64 = 5.0,    # energy cost of reproduction
    input_energy_density::Float64 = 100.0,
    # food properties
    spread_rate::Float64 = 0.9,
    spread_cost::Float64 = 10.0,
    σ_food_cap::Float64 = 5.0,
    σ_regen_rate::Float64 = 5.0,
)
    # sanity checks
    @assert size(food_distribution) == dims
    @assert lifetime > 0
    food_cap = maximum(food_distribution) * food_cap_multiplier

    # reproducibility matters
    rng = MersenneTwister(42)

    food_properties = (
        spread_rate = spread_rate,
        spread_cost = spread_cost,
        σ_food_cap = σ_food_cap,
        σ_regen_rate = σ_regen_rate,
    )

    # food lives in a separate model, since this improves performance
    food_space = GridSpace(dims; periodic = false)
    food_model = ABM(Food, food_space; properties = food_properties, rng)

    properties = (
        food = food_model,
        lifetime = lifetime,
        eat_rate_factor = eat_rate_factor,
        sensory_radius_cost = sensory_radius_cost,
        distance_cost = distance_cost,
        reproduction_energy_cost = reproduction_energy_cost,
        σ_sensory_radius = σ_sensory_radius,
        σ_reproduction_threshold = σ_reproduction_threshold,
        σ_speed = σ_speed,
        input_energy_density = input_energy_density,
    )

    # non-periodic euclidean space
    space = GridSpace(dims; periodic = false, metric = :euclidean)
    model = ABM(Bacterium, space; properties, rng)

    # add bacteria at random places
    for i ∈ 1:n_bacteria
        add_agent!(Bacterium, model, initial_bacterium)
    end

    # add food wherever it is nonzero
    for i = 1:dims[1], j = 1:dims[2]
        food_distribution[i, j] > 0 || continue

        add_agent!(
            (i, j),
            model.food,
            food_cap,
            food_distribution[i, j],
            initial_regen_rate,
        )
    end

    return model
end
