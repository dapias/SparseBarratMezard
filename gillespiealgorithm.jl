include("src/SparseBarratMezard.jl")
using JLD, JLD2, StatsBase
using .SparseBarratMezard
    
c = 5
T_grid = [0.2, 0.35, 0.5, 0.8]
#T_grid = [0.35, 0.5]
t_max = 10^4
ensemble = 10^4
for T in T_grid
    trees = Array{Dict{Float64,Node}, 1}(undef, ensemble);
    for i in 1:ensemble
        tree = init_tree(c, T)
        time= 0.0
        n_node, time = init_gillespie!(tree, c, T)
        n_node, time = gillespie1st!(tree, c, T, n_node, time);

        while time < t_max
            n_node, time = gillespie!(tree, c, T, n_node, time)
        end
        println(i)
        trees[i] = tree
    end
       
    @JLD2.save "data/gillespie/T=$(T)c=$(c).jld2" trees
end
