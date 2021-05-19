"""
    food_step!(model::ABM)

Simulates growth of food on the grid
"""
function food_step!(model::ABM)
    # randomise the update order for updating tiles
    shuffle!(model.rng, model.food_update_order)
    # update tiles
    for pos in model.food_update_order
        # if tile is empty
        if model.food[pos] == 0.0
            # regenerate food if the `random_spread_chance` occurs
            rand(model.rng) <= model.food_data.random_spread_chance &&
                (model.food[pos] += model.food_data.regen_rate; continue)

            # iterate through neighboring tiles
            for offset in model.food_data.neighborhood
                # bounds check
                all((0, 0) .< Tuple(pos + offset) .< size(model.food)) || continue
                # spread to this tile with probability proportional to amount of food
                # for neighboring dead tiles, this doesn't do anything
                if rand(model.rng) <=
                   model.food_data.spread_multiplier * model.food[pos+offset] /
                   model.food_data.food_cap
                    model.food[pos+offset] += model.food_data.regen_rate
                    break
                end
            end
        else
            # the tile is a food tile, so regenerate food till the cap
            model.food[pos] =
                min(model.food[pos] + model.food_data.regen_rate, model.food_data.food_cap)
        end
    end
end
