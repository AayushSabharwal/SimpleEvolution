"""
    reproduce!(bact::Bacterium, model::ABM)

Simulates reproduction of bacteria. Energy is halved after subtracting a cost. Inherited
features are varied.
"""
function reproduce!(bact::Bacterium, model::ABM)
    add_agent!(
        bact.pos,
        Bacterium,
        model,
        bact.species,
        0,
        (bact.energy - model.reproduction_energy_cost) / 2.0,
        inherit(bact.sensory_radius, model.rng, model.σ_sensory_radius),
        inherit(bact.reproduction_threshold, model.rng, model.σ_reproduction_threshold),
        inherit(bact.speed, model.rng, model.σ_speed),
        bact.food_target,
    )

    bact.age = 0
    bact.energy = (bact.energy - model.reproduction_energy_cost) / 2.0

    bact.sensory_radius = inherit(bact.sensory_radius, model.rng, model.σ_sensory_radius)
    bact.reproduction_threshold =
        inherit(bact.reproduction_threshold, model.rng, model.σ_reproduction_threshold)
    bact.speed = inherit(bact.speed, model.rng, model.σ_speed)
end

function agent_step!(bact::Bacterium, model::ABM)
    # die of age or starvation
    (bact.age > model.lifetime || bact.energy <= 0) && (kill_agent!(bact, model); return)

    # grow older
    bact.age += 1

    # reproduce if possible
    bact.energy >= bact.reproduction_threshold && reproduce!(bact, model)

    # if currently on a food source
    if model.food[bact.pos...] > 0.0
        # eat
        eaten = min(model.food[bact.pos...], bact.speed * model.eat_rate_factor)
        model.food[bact.pos...] -= eaten
        bact.energy += eaten
        return
    else
        # can't eat, so look for something else
        bact.food_target = (-1, -1)
    end

    # if bact doesn't see any food
    if bact.food_target == (-1, -1)
        best_pos = (-1, -1)
        # look for the best food nearby
        for pos in nearby_positions(bact.pos, model, bact.sensory_radius)
            best_pos == (-1, -1) || model.food[best_pos...] < model.food[pos...] || continue
            best_pos = pos
        end
        best_pos != (-1, -1) && (bact.food_target = best_pos)
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
        model.sensory_radius_cost * bact.sensory_radius +
        model.distance_cost * max(abs.(movement)...)
end
