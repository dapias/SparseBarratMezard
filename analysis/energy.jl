push!( LOAD_PATH, "/Users/diegotapias/cavity/src");
using SparseBarratMezard, JLD2, JLD, StatsBase

t_max = 10^4
t_grid = exp10.(range(0, stop = log10(t_max), length = 100));
T_grid = [0.2, 0.35, 0.5, 0.8]
#T_grid = [0.2, 0.35]
c = 5

for T in T_grid
    @JLD2.load "/Users/diegotapias/cavity/data/gillespie/trees/T=$(T)c=$(c).jld2" trees
    ensemble = length(trees)
    es = zeros(length(t_grid), ensemble);
    
    for i in 1:ensemble
        slices = Array{Node}(undef, length(t_grid));
        tt = sort(trees[i], by=x->x[1]);
        times = collect(keys(tt))
        nodes = collect(values(tt));
    for j in 1:length(t_grid)
        index = findlast(times.< t_grid[j])
        slices[j]  = nodes[index]
    end
    es[:,i] = [k.energy for k in slices]
    end

    meane = -[mean(es[k,:]) for k in 1:length(t_grid)]

     
    save("/Users/diegotapias/cavity/data/gillespie/energy/T=$(T)c=$(c).jld",  "t_grid", t_grid, "meane", meane, "ensemble", ensemble)
    println("T = $(T)")
end

