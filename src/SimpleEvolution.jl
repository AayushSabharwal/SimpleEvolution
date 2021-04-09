using DrWatson
using BSON

includet("core.jl")
includet("bacterium_step.jl")
includet("food_step.jl")
includet("runner.jl")
includet("plotter.jl")

food_distribution = reshape([max(600.0-(abs(i-20)*abs(j-20))^1.5, 0) for i in 1:40 for j in 1:40], 40, 40);
