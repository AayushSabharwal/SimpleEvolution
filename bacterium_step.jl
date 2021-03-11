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
    reproduce!(bact::Bacterium, model::ABM)

Simulates reproduction of bacteria. Energy is halved after subtracting a cost. Inherited features are
varied.
"""
function reproduce!(bact::Bacterium, model::ABM)
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
    bact.speed = inherit(bact.speed, model.rng, model.speed_std)
end

function agent_step!(bact::Bacterium, model::ABM)
    # die of age or starvation
    (bact.age > model.lifetime || bact.energy <= 0) && (kill_agent!(bact, model); return)

    # grow older
    bact.age += 1

    # reproduce if possible
    bact.energy >= bact.reproduction_threshold && reproduce!(bact, model)

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
