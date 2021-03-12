"""
    reproduce!(food::Food, pos::Dims{2}, model::ABM)

Simulates reproduction of food. Inherited features are varied.
"""
function reproduce!(food::Food, pos::Dims{2}, model::ABM)
    add_agent!(
        pos,
        Food,
        model,
        inherit(food.food_cap, model.rng, model.σ_food_cap),
        (food.current_food - model.spread_cost) / 2.0,
        inherit(food.regen_rate, model.rng, model.σ_regen_rate),
    )

    food.current_food = (food.current_food - model.spread_cost) / 2.0
end

@inline inbounds(pos::Dims{2}, dims::Dims{2}) = all(1.0 .<= pos .<= dims)

function agent_step!(food::Food, model::ABM)
    food.current_food <= 0 && (kill_agent!(food, model); return)

    food.current_food = min(food.current_food + food.regen_rate, food.food_cap)

    food.current_food > model.spread_cost || return

    for offset in NEIGHBORHOOD
        inbounds(food.pos .+ offset, size(model.space.s)) || continue
        isempty(food.pos .+ offset, model) || continue
        rand(model.rng) <= model.spread_rate || continue
        reproduce!(food, food.pos .+ offset, model)
        return
    end
end

function food_step!(model::ABM)
    Agents.step!(model.food, agent_step!)
end
