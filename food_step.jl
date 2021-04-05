function food_step!(model::ABM)
    shuffle!(model.rng, model.food_update_order)
    for pos in model.food_update_order
        if model.food[pos] == 0.0
            rand(model.rng) <= model.food_data.random_spread_chance && (model.food[pos] += model.food_data.regen_rate; continue)

            for offset in model.food_data.neighborhood
                all((0, 0) .< Tuple(pos+offset) .< size(model.food)) || continue
                if rand(model.rng) <= model.food_data.spread_multiplier * model.food[pos+offset]/model.food_data.food_cap
                    model.food[pos+offset] += model.food_data.regen_rate
                    break
                end
            end
        else
            model.food[pos] = min(model.food[pos] + model.food_data.regen_rate, model.food_data.food_cap)
        end
    end
end
