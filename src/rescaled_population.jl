export generate_populationr, population_updater!, DOSr
    
struct Population
    zetas::Array{Complex{Float64},1}
    energies::Array{Float64,1}
end

function symmetric_f(e1::Float64, beta::Float64, e2::Array{Float64})
    exp.(beta*(e1 .+ e2)/2)./(2*cosh.(beta*(e1.-e2)./2.)) 
end

function init_population(Np::Int64)
    dist_energy = Exponential()
    Population(rand(Complex{Float64},Np), rand(dist_energy, Np))
end

function generate_populationr(lambda::Float64,c::Int64,T::Float64,Np::Int64, nsteps::Int64, epsilon::Float64)
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
        zeta_cm1 += im*(lambda-epsilon*im)/(exp(-beta*e1)) 
        relement = rand(1:Np)
        poparray.zetas[relement] = zeta_cm1
        poparray.energies[relement] = e1
    end
    poparray
end

function population_updater!(lambda::Float64,c::Int64,T::Float64,Np::Int64, nsteps::Int64, epsilon::Float64, poparray::Population)
    beta = 1.0/T
    dist_energy = Exponential()
    for i in 1:nsteps
        random_elements = rand(1:Np,c-1)
        zetas = poparray.zetas[random_elements]
        energies = poparray.energies[random_elements]
        e1 = rand(dist_energy)
        fsym = symmetric_f(e1, beta, energies)
        zeta_cm1 = sum(im*fsym.*zetas./(im*fsym .+ zetas))
        zeta_cm1 += im*(lambda-epsilon*im)/(exp(-beta*e1)) 
        relement = rand(1:Np)
        poparray.zetas[relement] = zeta_cm1
        poparray.energies[relement] = e1
    end
    poparray
end

function DOSr(lambda::Float64, c::Int64, T::Float64, Np::Int64, ensemble::Int64, nsteps::Int64,  epsilon2::Array{Float64,1}; epsilon = 1.0e-300)
    poparray = init_population(Np)
    population_update!(lambda, c, T, Np, nsteps, epsilon, poparray);  ##Equilibrium is reached
    DOSr(lambda, c, T, Np, ensemble, epsilon2, poparray, epsilon = epsilon)
end

function DOSr(lambda::Float64, c::Int64, T::Float64, Np::Int64, ensemble::Int64, epsilon2::Array{Float64,1}, poparray::Population; epsilon = 1e-300)
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
        zeta_c =  im* (lambda-epsilon2[i]*im)/(exp(-beta*e1))  + sum_term
        res[1, i] = real(1/(zeta_c*(exp(-beta*e1))))
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
            zeta_c =   im*(lambda-epsilon2[i]*im)/(exp(-beta*e1))  + sum_term
            res[j, i] = real(1/(zeta_c*(exp(-beta*e1))))
        end
        ##update the population with new values
        population_update!(lambda,c,T,Np, 1, epsilon, poparray)
        ##########################
    end
    [mean(res[:,i])*1/pi for i in 1:length(epsilon2)]
end

function DOSr(l_grid::Union{Array, StepRangeLen}, c::Int64, T::Float64, Np::Int64, ensemble::Int64, nsteps::Int64,  epsilon2::Array{Float64,1})
    res = [DOSr(i, c, T, Np, ensemble, nsteps, epsilon2) for i in l_grid];
    res = vcat(res...);
end

