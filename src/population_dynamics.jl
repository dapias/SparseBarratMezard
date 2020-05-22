export generate_population, population_update!, DOS, DOSIPR, symmetric_f
    
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

function generate_population(lambda::Float64,c::Int64,T::Float64,Np::Int64, nsteps::Int64, epsilon::Float64)
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

function DOS(lambda::Float64, c::Int64, T::Float64, Np::Int64, ensemble::Int64, nsteps::Int64,  epsilon2::Array{Float64,1}; epsilon = 1.0e-300)
    poparray = init_population(Np)
    population_update!(lambda, c, T, Np, nsteps, epsilon, poparray);  ##Equilibrium is reached
    DOS(lambda, c, T, Np, ensemble, epsilon2, poparray, epsilon = epsilon)
end

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

