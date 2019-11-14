#using  Distributed
#addprocs(11)

@everywhere using JLD
@everywhere include("cavitymethod2.jl")
@everywhere function locincluster(lambda, T; N = 10^4, c = 3, epsilon2 =  [1e-5,1e-3, 1e-4])
    nsteps = N*10^4
    rho_mean = mean(DOS(lambda, c, T, 10^3, 10^5, 10^7, epsilon2))
    ns = 2^17
    #ns = 2^8
    epsi = 8.0./(rho_mean.*ns);

    Ms, energies, nei = generate_barrat_matrix(ns, c, 1/T)

    dos, ipr = DOSIPR(lambda, c, ns, 1/T, epsi, Ms, energies, nei)

   
    save("localization/scaling/single_instance/T=$(T)/lambda=$(lambda)N=$(N).jld", "n", ns, "Ms", Ms, "epsilon", epsi, "c", c, "energies", energies, "ipr", ipr)

    return ipr
end



@everywhere l_grid = collect(-1:0.01:0.)
@everywhere T_grid = 0.4:0.02:0.6
println("starting")

for T in T_grid
    println(T)
    #mkdir("localization/scaling/single_instance/T=$(T)/")
    ip = pmap(x -> locincluster(x,T), l_grid)
    save("localization/scaling/single_instance/ipr/T=$(T).jld", "ipr", ip, "lambda",l_grid)    
 end
