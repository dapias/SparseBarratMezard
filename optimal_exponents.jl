using JLD, Statistics

peak = ARGS[1]
orientation =ARGS[2]
c = ARGS[3]
#a = load("datos/collapse/density/$(peak)-$(orientation)-c=$(c)Bouchaud.jld")
a = load("datos/collapse/density/$(peak)-$(orientation)-c=$(c).jld")
#a = load("datos/collapse/density/$(peak)-$(orientation).jld")

peak = parse(Int64, peak)
rho = a["density"];
dlambda = a["parameters"]["dlambda"];
epsilons = a["parameters"]["eps2"];
Ts = a["parameters"]["T"];


xs_grid = collect(0.0:0.01:1.0);
xs_grid = repeat(xs_grid, 1, length(Ts))

dens_scaled = [epsilons[k].^xs_grid[:,1] for k in 1:length(epsilons)];
rho_scaled = [rho[:,1,:][:,k].*epsilons[k].^(xs_grid[:,1][1] ) for k in 1:4]
for i in 2:length(xs_grid[:,1])
    rho_temp = [rho[:,1,:][:,k].*epsilons[k].^(xs_grid[:,1][i] ) for k in 1:4]
    global rho_scaled = hcat(rho_scaled, rho_temp)
end

rho_diffTs = Array{Array{Float64,1},3}(undef, 4, length(xs_grid[:,1]), length(Ts));
rho_diffTs[:,:,1] = rho_scaled;
for zeta in 2:length(Ts)
    global rho_scaled = [rho[:,zeta,:][:,k].*epsilons[k].^(xs_grid[:,zeta][1] ) for k in 1:4]
    for i in 2:length(xs_grid[:,1])
        rho_temp = [rho[:,zeta,:][:,k].*epsilons[k].^(xs_grid[:,zeta][i] ) for k in 1:4]
        rho_scaled = hcat(rho_scaled, rho_temp)
    end
    global rho_diffTs[:,:,zeta] = rho_scaled;
end


means = [Statistics.mean([rho_diffTs[:,:,1][:,1][k][i] for k in 1:4]) for i in 1:31]
for exponente in 2:length(xs_grid[:,1])
    mean_temp =  [Statistics.mean([rho_diffTs[:,:,1][:,exponente][k][i] for k in 1:4]) for i in 1:31]
    global means = hcat(means, mean_temp)
end

rel_error = [sum([sum([rho_diffTs[:,:,1][:,exponente][k][i] .- means[i,exponente] for
                       k in 1:4].^2)./(means[i, exponente].^2) for i in 1:31]) for exponente in 1:length(xs_grid[:,1])];

for temperature in 2:length(Ts)
    means = [Statistics.mean([rho_diffTs[:,:,temperature][:,1][k][i] for k in 1:4]) for i in 1:31]
    for exponente in 2:length(xs_grid[:,1])
        mean_temp =  [Statistics.mean([rho_diffTs[:,:,temperature][:,exponente][k][i] for k in 1:4]) for i in 1:31]
        means = hcat(means, mean_temp)
    end
    
    temp_error = [sum([sum([rho_diffTs[:,:,temperature][:,exponente][k][i] .- means[i,exponente] for
                            k in 1:4].^2)./(means[i, exponente].^2) for i in 1:31]) for exponente in 1:length(xs_grid[:,1])];
    
    global rel_error = hcat(rel_error, temp_error)
end

ind_opt = [findmin(rel_error[:,i])[2] for i in 1:length(Ts)]
x_opt = [xs_grid[ind_opt[i], i] for i in 1:length(Ts)];

#save("datos/collapse/exponents/$(peak)-$(orientation)-c=$(c)Bouchaud.jld", "Ts", Ts, "xs", x_opt,   "rel_error", rel_error, "xs_grid", xs_grid)

save("datos/collapse/exponents/$(peak)-$(orientation)-c=$(c).jld", "Ts", Ts, "xs", x_opt,   "rel_error", rel_error, "xs_grid", xs_grid)

#save("datos/collapse/exponents/$(peak)-$(orientation).jld", "Ts", Ts, "xs", x_opt,
#     "rel_error", rel_error, "xs_grid", xs_grid)
