@everywhere parameters =
Dict(
"Np" => 1000,
"nsteps" => Int64(10^7),
"epsilon" => 1.e-300,
"c" => 5,
"ensemble" => Int64(10^6),
"eps2" => [1e-4, 1e-5, 1e-6, 1e-7],
    "T" => collect(0.05:0.1:0.95),
"dlambda" => exp10.(range( -1, stop = 2, length = 31)),
    "orientation" => "left")

peak = 0
population = generatepopulation(parameters, peak)
orientation = parameters["orientation"]
c = parameters["c"]
save("datos/collapse/population/$(peak)-$(orientation)-c=$(c).jld", "parameters", parameters, "peak", peak, "population", population)

@everywhere peak = 0
@everywhere orientation = "left"
@everywhere c = 5
@everywhere  a = load("datos/collapse/population/$(peak)-$(orientation)-c=$(c).jld")
@everywhere parameters = a["parameters"]
@everywhere population = a["population"]
@everywhere peak = a["peak"]


density = generatedensity(parameters, population, peak)

save("datos/collapse/density/$(peak)-$(orientation)-c=$(c).jld", "parameters", parameters, "peak", peak, "density", density)

#d2 = [density[:,1,:][:,k].*eps2[k]^(0.9) for k in 1:4]
