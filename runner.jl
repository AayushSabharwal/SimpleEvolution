using DataFrames
using Statistics

function run!(model::ABM; nsteps::Int = 2000)
    bact_data = DataFrame(step = Int[], id = Int[], energy = Float64[], sensory_radius = Float64[], repr_threshold = Float64[], speed = Float64[])
    food_data = DataFrame(step = Int[], id = Int[], x = Int[], y = Int[], food_cap = Float64[], current_food = Float64[], regen_rate = Float64[])

    for i in 0:nsteps
        for a in allagents(model)
            push!(bact_data, (i, a.id, a.energy, a.sensory_radius, a.reproduction_threshold, a.speed))
        end
        for a in allagents(model.food)
            push!(food_data, (i, a.id, a.pos[1], a.pos[2], a.food_cap, a.current_food, a.regen_rate))
        end
        step!(model, agent_step!, food_step!)
    end
    return bact_data, food_data
end
