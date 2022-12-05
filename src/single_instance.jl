export generate_sparse_barrat_matrix, generate_barrat_matrix, rho, DOSIPR, resolvent

###Generates the instances
function generate_sparse_barrat_elements(n::Int64,c::Int64)
    dene = Exponential()
    L = random_regular_graph(n,c)
    adj = adjacency_matrix(L);
    energies = rand(dene,nv(L));
    #M = 1.0*adj;
    l_nei = LightGraphs.SimpleGraphs.fadj(L)

    return energies, l_nei, adj, L
end

function generate_sparse_barrat_matrix(energies::Array, list_nei::Array, adj::SparseMatrixCSC, n::Int64, c::Int64, T::Float64)

    beta = 1/T
    M = 1.0*adj
    a = zeros(n)
    
    for i in 1:n
        for j in 1:c
            k = list_nei[i][j]
            M[i,k] = 1.0/c*(1.0/(1.0+exp(-beta*(energies[i] - energies[k]))))
            a[k] += M[i,k]
        end
    end
    
    for i in 1:n
        M[i,i] = -a[i]
    end
    
    peq = exp.(beta*energies)
    peq /= sum(peq)
    p_inv = peq.^(-1/2)
    p_dir = peq.^(1/2)
    Ms = Diagonal(p_inv)*M *Diagonal(p_dir);

end

function generate_barrat_matrix(energies::Array, list_nei::Array, adj::SparseMatrixCSC, n::Int64, c::Int64, T::Float64)

    Ms = generate_sparse_barrat_matrix(energies, list_nei, adj, n, c, T)
    Ms = Matrix(Symmetric(Ms))
    
    return Ms
end

####Cavity system
function omega_generator(adj::SparseMatrixCSC{Int64,Int64}, n::Int64, c::Int64)
    col = repeat(1:n, inner = c)
    sparse(adj.rowval, col, rand(Complex{Float64}, c*n), n, n)
end

function symmetric_f(e1::Float64, beta::Float64, e2::Float64)
    exp(beta*(e1 + e2)/2)/(2*cosh(beta*(e1-e2)/2.)) 
end

function local_error(lambda::Float64,  neighs::Array{Array{Int64,1},1},
        Omega_old::SparseMatrixCSC{Complex{Float64},Int64}, epsilon::Float64,  beta::Float64, energies::Array{Float64,1}, n::Int64, Omega_sum::SparseMatrixCSC{Complex{Float64},Int64}, c::Int64)

    for k in 1:n
        value = neighs[k]
        for j in value
            a = 0.0 +0.0*im
            for l in value
                if l != j
                    fs = symmetric_f(energies[k], beta, energies[l])
                    a+= im*fs*Omega_old[l,k]/(im*fs + Omega_old[l,k])
                end
            end
            @inbounds Omega_sum[k,j] = a + im*(lambda - im*epsilon)/(exp(-beta*energies[k])/c)
        end
    end
  
    error = sum(abs.(Omega_sum .- Omega_old))  

 #   @inbounds fill!(Omega_old, zero(Complex{Float64}))

    return error, Omega_sum, Omega_old
end

function cavity_omegas(lambda::Float64, c::Int64, epsilon::Float64,  energies::Array{Float64,1}, nei::Array{Array{Int64,1},1},
                       adj::SparseMatrixCSC{Int64,Int64}, T::Float64; tolerance = 0.1)
    ##The Barrat Matrix is not a needed  argument. The energies and the list of neighbours are enough information to apply the cavity method. The returned are the cavity_omegas

    beta = 1/T
    error2 = 10.0*tolerance*n*n
    Omegas = omega_generator(adj, n, c)
    Omegap = similar(Omegas)
    
    while error2 > tolerance
         error2, Omegas, Omegap = local_error(lambda, nei,  Omegas, epsilon,  beta, energies, n, Omegap, c)
    end
    
    return Omegas
end

function marginal_gs(lambda::Float64, c::Int64, epsilon::Float64, energies::Array{Float64,1}, nei::Array{Array{Int64,1},1}, n::Int64, Omegas::SparseMatrixCSC{Complex{Float64},Int64}, T::Float64)

    beta = 1/T
    mars = zeros(n)

    for j in 1:n
        sum_j = 0.0
        for k in nei[j]
            fs = symmetric_f(energies[k], beta, energies[j])
            @inbounds sum_j +=  im*fs*Omegas[k,j]/(im*fs + Omegas[k,j])
        end

        omega_j =  im*(lambda - im*epsilon)/(exp(-beta*energies[j])/c) + sum_j
        mars[j] = real(exp(beta*energies[j])/(c*omega_j))
    end

    mars
end

function dens(lambda::Float64, c::Int64, epsilon::Float64, energies::Array{Float64,1}, nei::Array{Array{Int64,1},1}, n::Int64, Omegas::SparseMatrixCSC{Complex{Float64},Int64}, T::Float64)

    sum(marginal_gs(lambda, c, epsilon, energies, nei, n, Omegas, T))/(pi*n)

end


        




function DOSIPR(lambda::Float64, c::Int64, n::Int64, beta::Float64, epsilon::Float64; tolerance = 0.1)
    A, energies, nei = generate_sparse_barrat_matrix(n,c, beta)
    DOSIPR(lambda, c, n, beta, epsilon, A, energies, nei)
end
    

function DOSIPR(lambda::Float64, c::Int64, beta::Float64, epsilon::Float64, energies::Array{Float64,1}, nei::Array{Array{Int64,1},1} ; tolerance = 0.1)
    n = length(energies)
    error2 = 10.0*tolerance*n*n
    Omegas = spzeros(Complex{Float64}, n, n)   ##Initialize
    fs = spzeros(Complex{Float64}, n, n)

    for key in 1:n
        value = nei[key]
        for i in value
            Omegas[key, i] = rand(Complex{Float64})
            fkeyi = symmetric_f(energies[i], beta, energies[key])
            fs[key,i] = fs[i,key] = fkeyi
        end
    end

    Omegap = spzeros(Complex{Float64}, n, n)

    
    while error2 > tolerance
        error2, Omegas, Omegap = local_error(lambda, nei,  Omegas, epsilon,  beta, energies, c, fs, n, Omegap)
    end
    
    sum_var = 0.
    sum_ipr = 0.
    for j in 1:n
        sum_j = 0.
        for k in nei[j]
            sum_j += im*Omegas[k,j]*fs[k,j]/( im*fs[k,j] + Omegas[k,j])
        end
        
        omega_j = im*(lambda - epsilon*im)/(exp(-beta*energies[j])/c)  + sum_j
        sum_ipr += abs2(1.0/(omega_j*(exp(-beta*energies[j])/c)))
        sum_var += real(1/(omega_j*(exp(-beta*energies[j])/c)))
    end
    
    dos = 1/(pi*n)*sum_var
    ipr = 1/n*sum_ipr
    
    return dos, epsilon*ipr/(dos*pi)
end





################# Analysis of the resolvent elements

function resolvent(lambda::Float64, c::Int64, beta::Float64, epsilon::Float64,  energies::Array{Float64,1}, nei::Array{Array{Int64,1},1}; tolerance = 0.1)
##The Barrat Matrix is not a needed  argument. The energies and the list of neighbours are enough information to apply the cavity method.
    n = length(energies)
    error2 = 10.0*tolerance*n*n
    Omegas = spzeros(Complex{Float64}, n, n)   
    fs = spzeros(Complex{Float64}, n, n)

    cavities = zeros(Float64, n)
    
    for key in 1:n
        value = nei[key]
        for i in value
            Omegas[key, i] = rand(Complex{Float64})
            fkeyi = symmetric_f(energies[i], beta, energies[key])
            fs[key,i] = fs[i,key] = fkeyi
        end
    end

    Omegap = spzeros(Complex{Float64}, n, n)

    while error2 > tolerance
         error2, Omegas, Omegap = local_error(lambda, nei,  Omegas, epsilon,  beta, energies, c, fs, n, Omegap)
    end
    
    sum_var = 0.
    for j in 1:n
        sum_j = 0.
        for k in nei[j]
            sum_j += im*Omegas[k,j]*fs[k,j]/( im*fs[k,j] + Omegas[k,j])
        end
        
        omega_j = im*(lambda - epsilon*im)/(exp(-beta*energies[j])/c) + sum_j
        #sum_var +=  real( (1/(omega_j*exp(-beta*energies[j])/c)))
        cavities[j] = real( (1/(omega_j*exp(-beta*energies[j])/c)))
    end

    return exp(mean(log.(cavities)))
    
end



