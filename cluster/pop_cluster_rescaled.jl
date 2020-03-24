#using  Distributed
#addprocs(3)

@everywhere include("../src/SparseBarratMezard.jl")
@everywhere using .SparseBarratMezard
@everywhere using JLD

@everywhere function pop_cluster(lambda, T; N = 10^4, c = 50, epsilon = 1e-300)
    nsteps = N*10^4
    pop_epsi =  generate_populationr(lambda, c, T, N, nsteps, epsilon)
    save("../data/populations_rescaled/T=$(T)/lambda=$(lambda)N=$(N)c=$(c)epsi_trescientos.jld", "n", N, "pop", pop_epsi)
    #save("../data/populations/T=$(T)/lambda=$(lambda)N=$(N)c=$(c)epsi_trescientos.jld", "n", N, "pop", pop_epsi)
end

@everywhere l_grid = collect(-3:0.02:0.)
#@everywhere l_grid = -exp10.(range(-3,stop =0, length = 101))
@everywhere T = 0.2

#mkdir("../data/populations/T=$(T)/")
println("start")
pmap(x -> pop_cluster(x,T), l_grid)

