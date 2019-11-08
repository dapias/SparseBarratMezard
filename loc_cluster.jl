using JLD
include("cavitymethod2.jl")
N = 10^2
nsteps = N*10^4
T = 0.6
lambda = -0.05
epsilon2 = [1e-5,1e-3, 1e-4];
rho_mean = mean(DOS(lambda, c, T, 10^3, 10^5, 10^7, epsilon2))
ns = 2^17
epsi = 8.0./(rho_mean.*ns);
pop_epsi =  population_update(lambda, c, T, N, nsteps, epsi)
save("localization/scaling/population/T=$(T)lambda=$(lambda)N=$(N)epsivariable.jld",
     "n", ns, "pop", pop_epsi)
