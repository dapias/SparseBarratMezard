export bouchaud_DOS, bouchaud_generate_population
    
struct Bouchaud_Population
    zetas::Array{Complex{Float64},1}
end

function bouchaud_init_population(Np::Int64)
    Bouchaud_Population(rand(Complex{Float64},Np))
end

function bouchaud_generate_population(lambda::Float64,c::Int64,T::Float64,Np::Int64, nsteps::Int64, epsilon::Float64)
    beta = 1/T
    poparray = bouchaud_init_population(Np)
    dist_energy = Exponential()
    for i in 1:nsteps
        random_elements = rand(1:Np,c-1)
        zetas = poparray.zetas[random_elements]
        #tau = rand(tau_dist)
        e1 = rand(dist_energy)
        zeta_cm1= im*(lambda-epsilon*im)/(exp(-beta*e1)/c)  + sum(im*zetas./(im .+ zetas))
        relement = rand(1:Np)
        poparray.zetas[relement] = zeta_cm1
    end
    poparray
end

function bouchaud_population_update!(lambda::Float64,c::Int64,T::Float64,Np::Int64, nsteps::Int64, epsilon::Float64, poparray::Bouchaud_Population)
    beta = 1.0/T
    dist_energy = Exponential()
    for i in 1:nsteps
        random_elements = rand(1:Np,c-1)
        zetas = poparray.zetas[random_elements]
        e1 = rand(dist_energy)
#        tau = rand(tau_dist)
        zeta_cm1= im*(lambda-epsilon*im)/(exp(-beta*e1)/c)  + sum(im*zetas./(im .+ zetas))
        relement = rand(1:Np)
        poparray.zetas[relement] = zeta_cm1
    end
    poparray
end

function bouchaud_DOS(lambda::Float64, c::Int64, T::Float64, Np::Int64, ensemble::Int64, nsteps::Int64,  epsilon2::Array{Float64,1}; epsilon = 1.0e-300)
    poparray = bouchaud_generate_population(lambda, c, T, Np,  nsteps, epsilon)
    bouchaud_DOS(lambda, c, T, Np, ensemble, epsilon2, poparray, epsilon = epsilon)
end

function bouchaud_DOS(lambda::Float64, c::Int64, T::Float64, Np::Int64, ensemble::Int64, epsilon2::Array{Float64,1}, poparray::Bouchaud_Population; epsilon = 1e-300)
    ##For this function a population array is passed that is assumed is already equilibrated
    beta = 1.0/T
    res = zeros(ensemble, length(epsilon2));
    dist_energy = Exponential()
    e1 = rand(dist_energy)
    random_elements = rand(1:Np,c)
    omega_sample = poparray.zetas[random_elements]
    for i  in 1:length(epsilon2)
        omega_c = im*(lambda-epsilon2[i]*im)/(exp(-beta*e1)/c) + sum(im*omega_sample./(1.0*im .+ omega_sample))
        res[1, i] = real(1/(omega_c*(exp(-beta*e1)/c)))
    end
    ##update the population after one measurement
    bouchaud_population_update!(lambda,c,T,Np, 1, epsilon, poparray)
    ##########################
    for j in 2:ensemble
        e1 = rand(dist_energy)
        random_elements = rand(1:Np,c)
        omega_sample = poparray.zetas[random_elements]
        for i  in 1:length(epsilon2)
            omega_c = im*(lambda-epsilon2[i]*im)/(exp(-beta*e1)/c) + sum(im*omega_sample./(1.0*im .+ omega_sample))
            res[j, i] = real(1/(omega_c*(exp(-beta*e1)/c)))
        end
        ##update the population with new values
        bouchaud_population_update!(lambda,c,T,Np, 1, epsilon, poparray)
        ##########################
    end
    [mean(res[:,i])*1/pi for i in 1:length(epsilon2)]
end

function bouchaud_DOS(l_grid::Union{Array, StepRangeLen}, c::Int64, T::Float64, Np::Int64, ensemble::Int64, nsteps::Int64,  epsilon2::Array{Float64,1})
    res = [bouchaud_DOS(i, c, T, Np, ensemble, nsteps, epsilon2) for i in l_grid];
    res = vcat(res...);
end

