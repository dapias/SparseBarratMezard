module SparseBarratMezard

using Distributions, LightGraphs, LinearAlgebra, SparseArrays, StatsBase

include("single_instance.jl")
include("population_dynamics.jl")
include("gillespie.jl")

end
