using Agents
using Random
using FileIO

const AXIAL_DIRECTIONS = [-1, 0, 1]
const NEIGHBORHOOD =
    [CartesianIndex(a) for a in Iterators.product([-1:1 for φ = 1:2]...) if a != (0, 0)]

@agent Bacterium GridAgent{2} begin
    species::Int
    age::Int        # age, in iterations. They only live so long
    energy::Float64 # energy. Less than zero, it dies. More than threshold, it doubles
    sensory_radius::Float64 # how far it can see food
    reproduction_threshold::Float64 # threshold energy to reproduce by binary fission
    speed::Int      # how far it can move in one iteration
    food_target::Dims{2}    # the food it is moving toward
end

struct FoodData
    neighborhood::Vector{CartesianIndex{2}} # what neighbors can food spread to
    food_cap::Float64   # limit on maximum food in one tile
    regen_rate::Float64 # regeneration rate for food tiles
    spread_multiplier::Float64  # proportionality constant for food spreading
    random_spread_chance::Float64   # probability of random spread
end

### Constructors
FoodData(;
    neighborhood::Vector{CartesianIndex{2}} = NEIGHBORHOOD,
    food_cap::Real = 200.0,
    regen_rate::Real = 10.0,
    spread_multiplier::Float64 = 0.07,
    random_spread_chance::Float64 = 0.0005,
) = FoodData(neighborhood, food_cap, regen_rate, spread_multiplier, random_spread_chance)

Bacterium(;
    energy::Real = 20.0,
    sensory_radius::Real = 4.0,
    reproduction_threshold::Real = 80.0,
    speed::Int = 3,
) = Bacterium(
    -1,
    (-1, -1),
    0,
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
    template.species,
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

"""
function initialize_model(food::Array{Float64,2}, n_bacteria::Union{Int,Vector{Int}}; kwargs...)

Initialize a model with given food distribution and number of bacteria

Keyword arguments
-----------------
- `initial_bacterium`: `Bacterium` denoting the initial values of parameters. `Vector{Bacterium}`
  for multiple species of bacteria. In such a case, `n_bacteria` should be a vector corresponding
  to the population of each species.
- `food_data` : `FoodData` parametrising the food spread
- `lifetime` : Maximum lifetime of a bacterium
- `eat_rate_factor` : Multiplied by a bacterium's `speed` to obtain the amount of food it will
  eat in one iteration
- `sensory_radius_cost`: Multipled by a bacterium's `sensory_radius` to obtain the corresponding
  energy cost
- `distance_cost` : Multipled by the distance a bacterium moves in one iteration to obtain the
  corresponding energy cost
- `σ_sensory_radius` : Standard deviation of normal distribution describing variation in
  `sensory_radius` over generations
- `σ_reproduction_threshold` : Similar to above, for `reproduction_threshold`
- `σ_speed` : Similar to above, for `speed`
- `reproduction_energy_cost` : Energy cost of reproduction
- `random_seed` : Seed for RNG used by the model
"""
function initialize_model(
    food::Array{Float64,2};
    n_bacteria::Union{Int,Vector{Int}} = 10,
    initial_bacterium::Union{Bacterium,Vector{Bacterium}} = Bacterium(),
    food_data::FoodData = FoodData(food_cap = maximum(food)),
    lifetime::Int = 200,
    eat_rate_factor::Float64 = 4.0,
    sensory_radius_cost::Float64 = 5.0,
    distance_cost::Float64 = 5.0,
    σ_sensory_radius::Float64 = 0.5,
    σ_reproduction_threshold::Float64 = 0.5,
    σ_speed::Float64 = 0.6,
    reproduction_energy_cost::Float64 = 10.0,
    random_seed::Int = 42,
)
    # sanity checks
    @assert lifetime > 0
    @assert (
        n_bacteria isa Array &&
        initial_bacterium isa Array &&
        length(n_bacteria) == length(initial_bacterium)
    ) || (!(n_bacteria isa Array) && !(initial_bacterium isa Array))

    dims = size(food)
    if initial_bacterium isa Array
        for i = 1:length(initial_bacterium)
            initial_bacterium[i].species = i
        end
    else
        initial_bacterium.species = 1
    end

    # reproducibility matters
    rng = MersenneTwister(random_seed)

    properties = (
        nspecies = initial_bacterium isa Array ? length(initial_bacterium) : 1,
        food = food,
        food_update_order = shuffle!(rng, collect(CartesianIndices(food))),
        food_data = food_data,
        lifetime = lifetime,
        eat_rate_factor = eat_rate_factor,
        sensory_radius_cost = sensory_radius_cost,
        distance_cost = distance_cost,
        reproduction_energy_cost = reproduction_energy_cost,
        σ_sensory_radius = σ_sensory_radius,
        σ_reproduction_threshold = σ_reproduction_threshold,
        σ_speed = σ_speed,
    )

    # non-periodic euclidean space
    space = GridSpace(dims; periodic = false, metric = :euclidean)
    model = ABM(Bacterium, space; properties, rng)

    # add bacteria at random places
    if n_bacteria isa Int
        for _ = 1:n_bacteria
            add_agent!(Bacterium, model, initial_bacterium)
        end
    else
        for species = 1:properties.nspecies
            for _ = 1:n_bacteria[species]
                add_agent!(Bacterium, model, initial_bacterium[species])
            end
        end
    end

    return model
end

"""
    image_to_foodmap(path::String)

Reads image at given path, returns a matrix representing relative food concentrations,
using the red channel of RGB value at each pixel.
"""
image_to_foodmap(path::String) = map(x -> x.r, load(path))
