####Gillespie

mutable struct Node  ##tengo que agregar la escape_rate
    distance::Int64
    visits::Int64
    energy::Float64
    energies::Array{Float64,1}
    parent::Union{Node, Nothing}
    rel_rates::Array{Float64,1}
    childs::Union{Array{Node,1}, Nothing}
    escape_rate::Float64
end

function node(init_node, i, beta)
    dene = Exponential()
    self_energy = init_node.energies[i]
    neighs_energy = vcat(init_node.energy, rand(dene,c-1))
    rates = 1.0/c.*(1.0./(1.0.+exp.(-beta.*(neighs_energy .- self_energy))));
    escape_rate = sum(rates)
    rel = rates./escape_rate;  
    Node(init_node.distance+1,0, self_energy, neighs_energy, init_node, rel, nothing, escape_rate)
end

function init_tree(c::Int64, T::Float64)
    dene = Exponential()   ##Energy distribution
    beta = 1/T
    self_energy = rand(dene)
    neighs_energy = rand(dene, c)

    rates = 1.0/c.*(1.0./(1.0.+exp.(-beta.*(neighs_energy .- self_energy))));
    escape_rate = sum(rates)
    #time_dist = Exponential(1.0./escape_rate)
    rel_rates = rates./escape_rate;
      
    initial_node = Node(0, 1,  self_energy, neighs_energy, nothing, rel_rates, nothing, escape_rate)

    time = 0.0
    tree = Dict(time => initial_node)
    
    return tree

end

function init_gillespie!(tree::Dict{Float64,Node}, c::Int64, T::Float64)
    beta = 1/T
    initial_node = tree[0.0]
    rel_rates = initial_node.rel_rates
    escape_rate = initial_node.escape_rate
    tw = rand(Exponential(1.0/escape_rate))

    ###crea los nodos vecinos
    childs = Array{Node}(undef, c)
    i = 1
    n_node = node(initial_node, i, beta)
    childs[1] = n_node
    #tree[n_node.distance] = [n_node]
    
    for i in 2:c
        n_node = node(initial_node, i, beta)
        childs[i] = n_node
        #ex = get!(tree, 1, [n_node])
        #push!(ex,n_node)
    end
    
    initial_node.childs = childs
 
    ####selecciona un nodo entre los vecinos
    next_node_index = sample(1:c, Weights(rel_rates));
    new_node = initial_node.childs[next_node_index]
    new_node.visits += 1
    tree[tw] = new_node

    return new_node, tw
end


function gillespie1st!(tree::Dict{Float64, Node},  c::Int64, T::Float64, current_node::Node,
                       time::Float64)
    
    dene = Exponential()
    beta = 1/T
    
    ###Los nuevos nodos tienen indice 2...c. El 1 es reservado al nodo padre
    childs = Array{Node}(undef, c-1)
    i = 2
    n_node = node(current_node, i, beta)
    childs[i - 1] = n_node
    for i in 3:c
        n_node = node(current_node, i, beta)
        childs[i-1] = n_node
        #ex = get!(tree, current_node.distance + 1, [n_node])
        #push!(ex,n_node)
    end
    current_node.childs = childs

    ###selecciona un nodo entre los vecinos
    escape_rate = current_node.escape_rate
    tw = rand(Exponential(1.0./escape_rate))
    rel_rates = current_node.rel_rates
    next_node_index = sample(1:c, Weights(rel_rates))

    if next_node_index == 1
        new_node = current_node.parent
    else
        new_node = current_node.childs[next_node_index-1]
    end
    new_node.visits += 1

    time += tw
    tree[time] = new_node

    return new_node, time
end


function gillespie!(tree::Dict{Float64, Node},  c::Int64, T::Float64, current_node::Node, time::Float64)
    dene = Exponential()
    beta = 1/T
    distance = current_node.distance
    rel_rates = current_node.rel_rates
    visits = current_node.visits
    
    escape_rate = current_node.escape_rate
    tw = rand(Exponential(1.0./escape_rate))
    next_node_index = sample(1:c, Weights(rel_rates))
    
    if current_node.childs !== nothing
        if distance == 0  ###Root node
            new_node = current_node.childs[next_node_index]
            new_node.visits += 1
        else
            if next_node_index == 1
                new_node = current_node.parent
                new_node.visits += 1
            else
                new_node = current_node.childs[next_node_index-1]
                new_node.visits += 1
            end
        end
    else   ###crea nuevos nodos
        childs = Array{Node}(undef, c-1)
        i = 2
        n_node = node(current_node, i, beta)
        childs[i - 1] = n_node
        for i in 3:c
            n_node = node(current_node, i, beta)
            childs[i-1] = n_node
        end
        current_node.childs = childs
        ###selecciona un nodo entre los vecinos
        if next_node_index == 1
            new_node = current_node.parent
        else
            new_node = current_node.childs[next_node_index-1]
        end
        new_node.visits += 1
    end

    time += tw
    tree[time] = new_node
    
    new_node, time
end

    
        
                


        

    
    

    


