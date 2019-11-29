using Distributions, LightGraphs, LinearAlgebra, SparseArrays


"""
    transition_rates(beta, c, energies, n::Int64)

For a given vector of energies gives a matrix whose entry (i,j) is equal to 1/c 1/(1+exp(-beta*(energies[i] - energies[j]))

"""
transition_rates(beta::Real, c::Int64, energies::Vector, n::Int64) = 1.0/c*[1.0/(1.0+exp(-beta*(energies[i] - energies[j]))) for i in 1:n, j in 1:n]

"""
    transition_rates(beta, c, energies, n::SimpleGraph)

For a given vector of energies gives a matrix whose entry (i,j) is equal to 1/c 1/(1+exp(-beta*(energies[i] - energies[j]))

"""
transition_rates(beta::Real, c::Int64, energies::Vector, grafo::SimpleGraph) = 1.0/c*[1.0/(1.0+exp(-beta*(energies[i] -
                                                                                                          energies[j]))) for i in vertices(grafo), j in vertices(grafo)]

function generate_sparse_barrat_matrix(n::Int64,c::Int64, beta::Float64)
    dene = Exponential()
    L = random_regular_graph(n,c)
    adj = adjacency_matrix(L);
    energies = rand(dene,nv(L));
    M = 1.0*adj;
    l_nei = L.fadjlist;

    for i in 1:n
        for j in 1:c
            k = l_nei[i][j]
            M[i,k] = 1.0/c*(1.0/(1.0+exp(-beta*(energies[i] - energies[k]))))
        end
    end

    for i in vertices(L)
        M[i,i] = -sum(M[:,i])
    end
    peq = exp.(beta*energies)
    peq /= sum(peq)
    p_inv = peq.^(-1/2)
    p_dir = peq.^(1/2)
    Ms = Diagonal(p_inv)*M *Diagonal(p_dir);


    return Ms, energies, l_nei
end

"""
    generate_barrat_matrix(n,c,beta)

Given the number of nodes (n) the connectivity (c) and the inverse of temperature (beta) generates a random matrix that corresponds to the symmetric version of the Barrat-Mézard model on a random regular graph with connectivity c. It also returns an array of the random energies associated with each node of the graph
"""
function generate_barrat_matrix(n::Int64,c::Int64, beta::Float64)

    Ms, energias, nei = generate_sparse_barrat_matrix(n,c,beta)

    Ms = Matrix(Symmetric(Ms))
    
    return Ms, energias, nei
end

#######################
####Single Instance

function symmetric_f(e1::Float64, beta::Float64, e2::Float64)
    exp(beta*(e1 + e2)/2)/(2*cosh(beta*(e1-e2)/2.)) 
end

function local_error(lambda::Float64,  neighs::Array{Array{Int64,1},1},
        Omega_old::SparseMatrixCSC{Complex{Float64},Int64}, epsilon::Float64,  beta::Float64, energies::Array, c::Int64, fs::SparseMatrixCSC{Complex{Float64},Int64}, N::Int64, Omega_sum::SparseMatrixCSC{Complex{Float64},Int64})

    reg = 10^-10.
    #N = length(energies)
    #Omega_sum = spzeros(Complex{Float64}, N, N)
    
    for k in 1:N
        value = neighs[k]
        for j in value
            for l in value
                if l != j
                    #fkl = symmetric_f(energies[l], beta, energies[k])
                    Omega_sum[k, j] += im*fs[k,l]*Omega_old[l, k]/(im*fs[k,l] +Omega_old[l,k])
                end
            end
            
            Omega_sum[k,j] += im*(lambda - im*epsilon)/(exp(-beta*energies[k])/c)
            
        end
    end

    relative_real = real.(Omega_sum .- Omega_old)./(abs.(real.(Omega_old))  .+ reg)
    relative_imag = imag.(Omega_sum .- Omega_old)./(abs.(imag.(Omega_old)) .+ reg)

    max_real = maximum(abs.(relative_real))
    max_imag = maximum(abs.(relative_imag))
    error = maximum([max_real, max_imag])
    #error = sum(abs.(Omega_sum .- Omega_old))

    @inbounds fill!(Omega_old, zero(Complex{Float64}))

    return error, Omega_sum, Omega_old
end

function rho(lambda::Float64, c::Int64, beta::Float64, epsilon::Float64, A::Union{Matrix, SparseMatrixCSC}, energies::
                        Array{Float64,1}, nei::Array{Array{Int64,1},1}; tolerance = 0.1)

    n = length(energies)
    error2 = 10.0*tolerance*n*n
    Omegas = spzeros(Complex{Float64}, n, n)   
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
    for j in 1:n
        sum_j = 0.
        for k in nei[j]
            #fkj = symmetric_f(energies[j], beta, energies[k])
            sum_j += im*Omegas[k,j]*fs[k,j]/( im*fs[k,j] + Omegas[k,j])
        end
        
        omega_j = im*(lambda - epsilon*im)/(exp(-beta*energies[j])/c) + sum_j
        sum_var +=  real( (1/(omega_j*exp(-beta*energies[j])/c)))
    end
    
    return 1/(pi*n)*sum_var
    
end

function rho(lambda::Float64, c::Int64, n::Int64, beta::Float64, epsilon::Float64; tolerance = 0.1)
    A, energies, nei  = generate_sparse_barrat_matrix(n,c, beta)

    rho(lambda, c,  beta, epsilon, A, energies, nei)
    
end


function DOSIPR(lambda::Float64, c::Int64, n::Int64, beta::Float64, epsilon::Float64; tolerance = 0.1)
    A, energies, nei = generate_sparse_barrat_matrix(n,c, beta)

    DOSIPR(lambda, c, n, beta, epsilon, A, energies, nei)
end
    

function DOSIPR(lambda::Float64, c::Int64, n::Int64, beta::Float64, epsilon::Float64,
                A, energies, nei; tolerance = 0.1)
    
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
            #fkj = symmetric_f(energies[j], beta, energies[k])
            sum_j += im*Omegas[k,j]*fs[k,j]/( im*fs[k,j] + Omegas[k,j])
        end
        
        omega_j = im*(lambda - epsilon*im)/(exp(-beta*energies[j])/c)  + sum_j
        sum_ipr += norm(1.0/(omega_j*(exp(-beta*energies[j])/c)))^2
        sum_var += real(1/(omega_j*(exp(-beta*energies[j])/c)))
    end
    
    dos = 1/(pi*n)*sum_var
    ipr = 1/n*sum_ipr
    
    return dos, epsilon*ipr/(dos*pi)
end


function omegas(lambda::Float64, c::Int64, beta::Float64, epsilon::Float64, A::Union{Matrix, SparseMatrixCSC}, energies::
                        Array{Float64,1}, nei::Array{Array{Int64,1},1}; tolerance = 0.1)

    n = length(energies)
    error2 = 10.0*tolerance*n*n
    Omegas = spzeros(Complex{Float64}, n, n)   ##Initialize
    
    for key in 1:n
        value = nei[key]
        for i in value
            Omegas[key, i] = rand(Complex{Float64})
        end
    end
    
    while error2 > tolerance
        error2, Omegas = local_error(lambda, nei,  Omegas, epsilon,  beta, energies, c)
    end
    
    oms = zeros(Complex{Float64}, n)
    for j in 1:n
        sum_j = 0.
        for k in nei[j]
            fkj = symmetric_f(energies[j], beta, energies[k])
            sum_j += im*Omegas[k,j]*fkj/( im*fkj + Omegas[k,j])
        end
        
        omega_j = im*(lambda - epsilon*im)/(exp(-beta*energies[j])/c) + sum_j
        oms[j] = 1/(omega_j*(exp(-beta*energies[j])/c))
    end
    
    return oms
    
end


#########################################
####Population Dynamics

"""
    Population

This struct is useful to represent any element of the ensemble evolved with the Population Dynamics algorithm. It consists of the pair (z, E)
"""
struct Population
    zetas::Array{Complex{Float64},1}
    energies::Array{Float64,1}
end



function symmetric_f(e1::Float64, beta::Float64, e2::Array{Float64})
    exp.(beta*(e1 .+ e2)/2)./(2*cosh.(beta*(e1.-e2)./2.)) 
end


"""
    init_population(Np::Int64)

Initializes a population {(z, E)} of size Np.
"""
function init_population(Np::Int64)
    dist_energy = Exponential()
    Population(rand(Complex{Float64},Np), rand(dist_energy, Np))
end

"""
    population_update(λ, c, T, Np, nsteps, ϵ)

Create and update the population {(z, E)} according to the cavity equations of the Barrat-Mézard master operator evaluated at λ for the set of parameters passed. 
"""
function population_update(lambda::Float64,c::Int64,T::Float64,Np::Int64, nsteps::Int64, epsilon::Float64)
    beta = 1.0/T
    poparray = init_population(Np)
    dist_energy = Exponential()
    for i in 1:nsteps
        random_elements = rand(1:Np,c-1)
        zetas = poparray.zetas[random_elements]
        energies = poparray.energies[random_elements]
        e1 = rand(dist_energy)
        fsym = symmetric_f(e1, beta, energies)
        zeta_cm1 = sum(im*fsym.*zetas./(im*fsym .+ zetas))
        zeta_cm1 += im*(lambda-epsilon*im)/(exp(-beta*e1)/c) 
        relement = rand(1:Np)
        poparray.zetas[relement] = zeta_cm1
        poparray.energies[relement] = e1
    end
    poparray
end

"""
    population_update!(λ, c, T, Np, nsteps, ϵ, poparray)

Update the population passed (set of {(ω, E)}) for the Barrat-Mézard master operator evaluated at λ for the set of parameters passed. 
"""
function population_update!(lambda::Float64,c::Int64,T::Float64,Np::Int64, nsteps::Int64, epsilon::Float64, poparray::Population)
    beta = 1.0/T
    dist_energy = Exponential()
    for i in 1:nsteps
        random_elements = rand(1:Np,c-1)
        zetas = poparray.zetas[random_elements]
        energies = poparray.energies[random_elements]
        e1 = rand(dist_energy)
        fsym = symmetric_f(e1, beta, energies)
        zeta_cm1 = sum(im*fsym.*zetas./(im*fsym .+ zetas))
        zeta_cm1 += im*(lambda-epsilon*im)/(exp(-beta*e1)/c) 
        relement = rand(1:Np)
        poparray.zetas[relement] = zeta_cm1
        poparray.energies[relement] = e1
    end
    poparray
end

"""
    DOS(lambda, c, T, Np, ensemble, nsteps, epsilon2::Array{Float64,1})

Create a population of size Np, update its elements according to the cavity equations for a number of steps equal to nsteps. Then extract elements of the population to compute the average <e^(β E) c/Ω_c> using a number of "ensemble" realizations. The value of epsilon is used to update the population up to equilibrium whereas the value of epsilon2 is used to compute the spectral density with the Lorentzians of width epsilon2.

"""
function DOS(lambda::Float64, c::Int64, T::Float64, Np::Int64, ensemble::Int64, nsteps::Int64,  epsilon2::Array{Float64,1}; epsilon = 1.0e-300)
    beta = 1.0/T
    poparray = init_population(Np)
    population_update!(lambda, c, T, Np, nsteps, epsilon, poparray);  ##Equilibrium is reached
    res = zeros(ensemble, length(epsilon2))
    dist_energy = Exponential()
    random_elements = rand(1:Np,c)
    zetas_sample = poparray.zetas[random_elements]
    energies = poparray.energies[random_elements]  
    e1 = rand(dist_energy)
    fsym = symmetric_f(e1, beta, energies)
    sum_term = sum(im*fsym.*zetas_sample./(im*fsym .+ zetas_sample))
    for i  in 1:length(epsilon2)
        zeta_c =   im*(lambda-epsilon2[i]*im)/(exp(-beta*e1)/c)  + sum_term
        res[1, i] = real(1/(zeta_c*(exp(-beta*e1)/c)))
    end
    ##update the population after one measurement
    population_update!(lambda,c,T,Np, 1, epsilon, poparray)
    ##########################
    for j in 2:ensemble
        random_elements = rand(1:Np,c)
        zetas_sample = poparray.zetas[random_elements]
        energies = poparray.energies[random_elements]  
        e1 = rand(dist_energy)
        fsym = symmetric_f(e1, beta, energies)
        sum_term = sum(im*fsym.*zetas_sample./(im*fsym .+ zetas_sample))
        for i  in 1:length(epsilon2)
            zeta_c =   im*(lambda-epsilon2[i]*im)/(exp(-beta*e1)/c)  + sum_term
            res[j, i] = real(1/(zeta_c*(exp(-beta*e1)/c)))
        end
        ##update the population with new values
        population_update!(lambda,c,T,Np, 1, epsilon, poparray)
        ##########################
    end
    [mean(res[:,i])*1/pi for i in 1:length(epsilon2)]
    
end


"""
    DOS(lambda, c, T, Np, ensemble, epsilon2, poparray::Population)

This function takes as the last argument a population that is assumed to be in equilibrium, then extract its elements to compute the average <e^(β E) c/Ω_c> using a number of "ensemble" realizations. The value of epsilon is used to update the population after a sample of e^(βE)c/Ω_c is taken whereas the value of epsilon2 is used to compute the spectral density with the Lorentzians of width epsilon2.
"""
function DOS(lambda::Float64, c::Int64, T::Float64, Np::Int64, ensemble::Int64, epsilon2::Array{Float64,1}, poparray::Population; epsilon = 1e-300)
    ##For this function a population array is passed that is assumed is already equilibrated

    beta = 1.0/T
    res = zeros(ensemble, length(epsilon2));
    dist_energy = Exponential()
    random_elements = rand(1:Np,c)
    zetas_sample = poparray.zetas[random_elements]
    energies = poparray.energies[random_elements]  
    e1 = rand(dist_energy)
    fsym = symmetric_f(e1, beta, energies)
    sum_term = sum(im*fsym.*zetas_sample./(im*fsym .+ zetas_sample))
    for i  in 1:length(epsilon2)
        zeta_c =  im* (lambda-epsilon2[i]*im)/(exp(-beta*e1)/c)  + sum_term
        res[1, i] = real(1/(zeta_c*(exp(-beta*e1)/c)))
    end
    ##update the population after one measurement
    population_update!(lambda,c,T,Np, 1, epsilon, poparray)
    ##########################
    for j in 2:ensemble
        random_elements = rand(1:Np,c)
        zetas_sample = poparray.zetas[random_elements]
        energies = poparray.energies[random_elements]  
        e1 = rand(dist_energy)
        fsym = symmetric_f(e1, beta, energies)
        sum_term = sum(im*fsym.*zetas_sample./(im*fsym .+ zetas_sample))
        for i  in 1:length(epsilon2)
            zeta_c =   im*(lambda-epsilon2[i]*im)/(exp(-beta*e1)/c)  + sum_term
            res[j, i] = real(1/(zeta_c*(exp(-beta*e1)/c)))
        end
        ##update the population with new values
        population_update!(lambda,c,T,Np, 1, epsilon, poparray)
        ##########################
    end
    [mean(res[:,i])*1/pi for i in 1:length(epsilon2)]
end

function DOS(l_grid::Union{Array, StepRangeLen}, c::Int64, T::Float64, Np::Int64, ensemble::Int64, nsteps::Int64,  epsilon2::Array{Float64,1})
    res = [DOS(i, c, T, Np, ensemble, nsteps, epsilon2) for i in l_grid];
    res = vcat(res...);
end



#######################################3
##################################
#################################


###Paso la poblacion equilibrada y computo simultaneamente DOS y IPR
function DOSIPR(lambda::Float64, c::Int64, T::Float64, Np::Int64, ensemble::Int64, epsilon2::Array{Float64,1}, poparray::Population; epsilon = 1e-300)
    ##For this function a population array is passed that is assumed is already equilibrated

    beta = 1.0/T
    ##init arrays
    res = zeros(ensemble, length(epsilon2));
    ipr = zeros(ensemble, length(epsilon2))
    ######
    dist_energy = Exponential()
    random_elements = rand(1:Np,c)
    zetas_sample = poparray.zetas[random_elements]
    energies = poparray.energies[random_elements]  
    e1 = rand(dist_energy)
    fsym = symmetric_f(e1, beta, energies)
    sum_term = sum(im*fsym.*zetas_sample./(im*fsym .+ zetas_sample))
    for i  in 1:length(epsilon2)
        zeta_c =   im*(lambda-epsilon2[i]*im)/(exp(-beta*e1)/c)  + sum_term
        res[1, i] = real(1/(zeta_c*(exp(-beta*e1)/c)))
        ipr[1,i] = (norm(1/(zeta_c*(exp(-beta*e1)/c))))^2
    end
    ##update the population after one measurement
    population_update!(lambda,c,T,Np, 1, epsilon, poparray)
    ##########################
    for j in 2:ensemble
        random_elements = rand(1:Np,c)
        zetas_sample = poparray.zetas[random_elements]
        energies = poparray.energies[random_elements]  
        e1 = rand(dist_energy)
        fsym = symmetric_f(e1, beta, energies)
        sum_term = sum(im*fsym.*zetas_sample./(im*fsym .+ zetas_sample))
        for i  in 1:length(epsilon2)
            zeta_c =   im*(lambda-epsilon2[i]*im)/(exp(-beta*e1)/c)  + sum_term
            res[j, i] = real(1/(zeta_c*(exp(-beta*e1)/c)))
            ipr[j,i] = (norm(1/(zeta_c*(exp(-beta*e1)/c))))^2
        end
        ##update the population with new values
        population_update!(lambda,c,T,Np, 1, epsilon, poparray)
        ##########################
    end
    dos = [mean(res[:,i])*1/pi for i in 1:length(epsilon2)]
    i2 = [mean(ipr[:,i])*epsilon2[i] for i in 1:length(epsilon2)]./(dos*pi)
    
    dos, i2
end

###Sin pasar la poblacion
function DOSIPR(lambda::Float64, c::Int64, T::Float64, Np::Int64, ensemble::Int64, nsteps::Int64, epsilon2::Array{Float64,1}; epsilon = 1e-300)

    beta = 1.0/T
    ###
    poparray = init_population(Np)
    population_update!(lambda, c, T, Np, nsteps, epsilon, poparray); 
    ##init arrays
    res = zeros(ensemble, length(epsilon2));
    ipr = zeros(ensemble, length(epsilon2))
    ######
    dist_energy = Exponential()
    random_elements = rand(1:Np,c)
    zetas_sample = poparray.zetas[random_elements]
    energies = poparray.energies[random_elements]  
    e1 = rand(dist_energy)
    fsym = symmetric_f(e1, beta, energies)
    sum_term = sum(im*fsym.*zetas_sample./(im*fsym .+ zetas_sample))
    for i  in 1:length(epsilon2)
        zeta_c =   im*(lambda-epsilon2[i]*im)/(exp(-beta*e1)/c)  + sum_term
        res[1, i] = real(1/(zeta_c*(exp(-beta*e1)/c)))
        ipr[1,i] = (norm(1/(zeta_c*(exp(-beta*e1)/c))))^2
    end
    ##update the population after one measurement
    population_update!(lambda,c,T,Np, 1, epsilon, poparray)
    ##########################
    for j in 2:ensemble
        random_elements = rand(1:Np,c)
        zetas_sample = poparray.zetas[random_elements]
        energies = poparray.energies[random_elements]  
        e1 = rand(dist_energy)
        fsym = symmetric_f(e1, beta, energies)
        sum_term = sum(im*fsym.*zetas_sample./(im*fsym .+ zetas_sample))
        for i  in 1:length(epsilon2)
            zeta_c =   im*(lambda-epsilon2[i]*im)/(exp(-beta*e1)/c)  + sum_term
            res[j, i] = real(1/(zeta_c*(exp(-beta*e1)/c)))
            ipr[j,i] = (norm(1/(zeta_c*(exp(-beta*e1)/c))))^2
        end
        ##update the population with new values
        population_update!(lambda,c,T,Np, 1, epsilon, poparray)
        ##########################
    end
    dos = [mean(res[:,i])*1/pi for i in 1:length(epsilon2)]
    i2 = [mean(ipr[:,i])*epsilon2[i] for i in 1:length(epsilon2)]./(dos*pi)
    
    dos, i2
end

