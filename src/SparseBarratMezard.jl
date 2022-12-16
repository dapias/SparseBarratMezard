module SparseBarratMezard

using Distributions, Graphs, LinearAlgebra, SparseArrays, StatsBase,
    InverseLaplace

include("single_instance.jl")
include("population_dynamics.jl")
#include("rescaled_population.jl")
#include("laplacetransform.jl")

end
