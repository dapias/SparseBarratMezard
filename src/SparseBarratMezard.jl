module SparseBarratMezard

using Distributions, LightGraphs, LinearAlgebra, SparseArrays, StatsBase,
    InverseLaplace

include("single_instance.jl")
#include("rescaled_single_instance.jl")
#include("bouchaud_single_instance.jl")
include("population_dynamics.jl")
include("rescaled_population.jl")
include("bouchaud_population_dynamics.jl")
include("gillespie.jl")
#include("bouchaud_gillespie.jl")
include("laplacetransform.jl")

end
