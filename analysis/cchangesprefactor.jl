@everywhere include("../src/SparseBarratMezard.jl")
@everywhere using .SparseBarratMezard, StatsBase
@everywhere using JLD, JLD2

@everywhere function cdependence(c, T; Np = 10^4, nsteps = Np*10^4, ensemble = 10^7, epsilon2 = [1e-5])
    
    lambda_barrat = -10^-2.5
    rho = [DOS(lambda_barrat, c, T, Np, ensemble, nsteps, epsilon2)[1] for i in 1:3]

 end

   
@everywhere connectivity = round.(Int64, exp10.(range(log10(5), stop = 2, length=10)))
@everywhere T = 0.35

dos = pmap(x -> cdependence(x, T), connectivity)

save("../data/dos/cdependenceT=$(T)mod.jld", "dos", dos, "cs", connectivity, "lambda",  -10^-2.5, "epsilon", 1e-5)
