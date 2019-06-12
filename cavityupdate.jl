#include("cavitymethod.jl")
#include("localization.jl")
using Distributions, JLD

struct Population2
    zetas::Array{Complex{Float64},1}
    energies::Array{Float64,1}
end

function barrat_rates(e_1, beta, es)
    1.0./(1.0 .+ exp.(-beta.*(es .- e_1)))
end

function asym_barrat_rates(e_1, beta, es)
    1.0 - big(1/(1 + exp(-beta*(e_1 - es))))
end

function init_population(Np)
    dist_energy = Exponential()
    Population2(rand(Complex{Float64},Np), rand(dist_energy, Np))
end


function pop_update(lambda,c,T,Np, nsteps, epsilon)
    beta = 1.0/T
    poparray = init_population(Np)
    dist_energy = Exponential()
    for i in 1:nsteps
        random_elements = rand(1:Np,c-1)
        zetas = poparray.zetas[random_elements]
        energies = poparray.energies[random_elements]
        e1 = rand(dist_energy)
        #Wj1 = big.(barrat_rates(e1, beta, energies))
        Wj1 = barrat_rates(e1, beta, energies)
        # for i in 1:length(Wj1)
        #     if Wj1[i] > 0.99
        #         Wj1[i] = asym_barrat_rates(e1, beta, energies[i])
        #     end
        # end
        W1j =  (1.0 .- Wj1)
        zeta_cm1 =  im*(lambda-epsilon*im)*c + sum(im*(Wj1.* zetas)./(im*(W1j) .+ zetas))
        relement = rand(1:Np)
        poparray.zetas[relement] = zeta_cm1
        poparray.energies[relement] = e1
    end
    poparray
end

function pop_update!(lambda,c,T,Np, nsteps, epsilon, poparray)
    beta = 1.0/T
    dist_energy = Exponential()
    for i in 1:nsteps
        random_elements = rand(1:Np,c-1)
        zetas = poparray.zetas[random_elements]
        energies = poparray.energies[random_elements]
        e1 = rand(dist_energy)
        #Wj1 = big.(barrat_rates(e1, beta, energies))
        Wj1 = barrat_rates(e1, beta, energies)
        # for i in 1:length(Wj1)
        #     if Wj1[i] > 0.99
        #         Wj1[i] = asym_barrat_rates(e1, beta, energies[i])
        #     end
        # end
        W1j =  (1.0 .- Wj1)
        zeta_cm1 =  im*(lambda-epsilon*im)*c + sum(im*(Wj1.* zetas)./(im*(W1j) .+ zetas))
        relement = rand(1:Np)
        poparray.zetas[relement] = zeta_cm1
        poparray.energies[relement] = e1
    end
    poparray
end


function sampling_equilibrium(lambda, c, T, Np, ensemble, epsilon, poparray)
    beta = 1.0/T
    new_population = zeros(Complex{Float64}, ensemble + Np)
    new_population[end-Np+1:end] = poparray.zetas
    dist_energy = Exponential()
    #Wj1 = zeros(BigFloat, c-1)
    for j in 1:ensemble
        random_elements = rand(1:Np,c-1)
        zetas_sample = poparray.zetas[random_elements]
        energies = poparray.energies[random_elements]
        e1 = rand(dist_energy)
        ##update the population with new values
        Wj1 = barrat_rates(e1, beta, energies)
        # for i in 1:length(Wj1)
        #     Wj1[i] = barrat_rates(e1, beta, energies[i])
        #     if Wj1[i] > 0.99
        #         Wj1[i] = asym_barrat_rates(e1, beta, energies[i])
        #     end
        # end
        W1j = 1.0 .- Wj1
        zeta_cm1 = im*(lambda-epsilon*im)*c +sum(im*(Wj1.* zetas_sample)./(im*(W1j) .+ zetas_sample))
        if real(zeta_cm1) > 10^50.
            println("Wj1", Wj1, "W1j", W1j, "zetas", zetas_sample,
                    "zeta_cm1", zeta_cm1)
        end
        relement = rand(1:Np)
        poparray.zetas[relement] = zeta_cm1
        poparray.energies[relement] = e1
        new_population[j] = zeta_cm1
        ##########################
    end
    new_population
    #[mean(res[:,i])*1/pi for i in 1:length(epsilons)], new_elements
end
    

function rho_barrat(lambda, c, T, Np, ensemble, epsilon, epsilons::Array{Float64,1}, poparray::Population2)
    ##For this function a population array is passed that is assumed is already equilibrated
    beta = 1.0/T
    marginal_variances_real = zeros(ensemble, length(epsilons));
    marginal_precisions = zeros(Complex{Float64}, ensemble, length(epsilons))
     new_population = zeros(Complex{Float64}, ensemble + Np)
    new_population[end-Np+1:end] = poparray.zetas
    dist_energy = Exponential()
    random_elements = rand(1:Np,c)
    zetas_sample = poparray.zetas[random_elements]
    energies = poparray.energies[random_elements]  
    e1 = rand(dist_energy)
    Wj1 = barrat_rates(e1, beta, energies)
    W1j = (1.0 .- Wj1)
    sum_term = sum(im*(Wj1.* zetas_sample)./(im*(W1j) .+ zetas_sample))
    for i  in 1:length(epsilons)
        zeta_c =  im*(lambda-epsilons[i]*im)*c + sum_term
        marginal_variances_real[1, i] = real(c/zeta_c)
        marginal_precisions[1, i] = zeta_c
    end
    ##update the population with new values
    zetas_sample = zetas_sample[1:end-1]
    energies = energies[1:end-1]
    Wj1 = barrat_rates(e1, beta, energies)
    taus = exp.(beta.*energies)
    W1j = 1.0 .- Wj1
    zeta_cm1 = im*(lambda-epsilon*im)*c +sum(im*(Wj1.* zetas_sample)./(im*(W1j) .+ zetas_sample))
    relement = rand(1:Np)
    poparray.zetas[relement] = zeta_cm1
    poparray.energies[relement] = e1
    new_population[1] = zeta_cm1
    ###################
    
    for j in 2:ensemble
        random_elements = rand(1:Np,c)
        zetas_sample = poparray.zetas[random_elements]
        energies = poparray.energies[random_elements]  
        e1 = rand(dist_energy)
        Wj1 = barrat_rates(e1, beta, energies)
        W1j =  (1.0 .- Wj1)
        sum_term = sum(im*(Wj1.* zetas_sample)./(im*(W1j) .+ zetas_sample))
        for i  in 1:length(epsilons)
            zeta_c =  im*(lambda-epsilons[i]*im)*c + sum_term
            marginal_variances_real[j, i] = real(c/zeta_c)
            marginal_precisions[j,i] = zeta_c
        end
        ##update the population with new values
        zetas_sample = zetas_sample[1:end-1]
        energies = energies[1:end-1]
        #Wj1 = big.(barrat_rates(e1, beta, energies))
        Wj1 = barrat_rates(e1, beta, energies)
        W1j = 1.0 .- Wj1
        zeta_cm1 = im*(lambda-epsilon*im)*c +sum(im*(Wj1.* zetas_sample)./(im*(W1j) .+ zetas_sample))
        relement = rand(1:Np)
        poparray.zetas[relement] = zeta_cm1
        poparray.energies[relement] = e1
        new_population[j] = zeta_cm1
        ##########################
    end
    
    marginal_variances_real./pi, new_population, marginal_precisions
end



