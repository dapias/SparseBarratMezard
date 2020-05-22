#using  Distributed
#addprocs(3)

@everywhere include("src/SparseBarratMezard.jl")
@everywhere using .SparseBarratMezard
@everywhere using JLD, StatsBase
@everywhere function locincluster(lambda, T; N = 10^4, c = 3, epsilon2 =  [1e-5,1e-3, 1e-4])
    nsteps = N*10^4
    rho_mean = mean(DOS(lambda, c, T, 10^3, 10^5, 10^7, epsilon2))
    n = [2^16, 2^17, 2^18]
    for ns in n
        epsi = 8.0./(rho_mean.*ns);
        pop_epsi =  generate_population(lambda, c, T, N, nsteps, epsi)
        save("localization/scaling/population/T=$(T)lambda=$(lambda)N=$(N)n=$(ns).jld",
             "n", ns, "pop", pop_epsi)
    end
end


@everywhere function pop_cluster(lambda, T; N = 10, c = 5, epsilon = 1e-300)
    nsteps = N*10^4
    pop_epsi =  generate_population(lambda, c, T, N, nsteps, epsi)
    save("data/population/T=$(T)/lambda=$(lambda)N=$(N)c=$(c)epsi_trescientos.jld",
         "n", N, "pop", pop_epsi)
end

@everywhere l_grid = collect(-0.4:0.01:-0.1)
@everywhere T = 0.2

#mkdir("data/population/T=$(T)/")
#pmap(x -> pop_cluster(x,T), l_grid)
pmap(x->locincluster(x, T), l_grid)
