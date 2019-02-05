include("cavitymethod.jl")
using PyCall

@pyimport scipy.special as sp

function meanK(ene, beta)
    exp(ene*beta)*sp.hyp2f1(1.0,1/beta, 1+1/beta, -exp(ene*beta))
end

function rho_barrat_cinfty(lambda, c, T, ensemble, epsilon, epsilon2)
    beta = 1.0/T
    res = zeros(ensemble);
    dist_energy = Exponential()
    obar = zeros(Complex{Float64},c)
    el = zeros(c)
    for i in 1:c
        e1 = rand(dist_energy)
        el[i] = e1
        K = meanK(e1, beta)
        obar[i] =  im*(lambda-epsilon*im)*exp(beta*e1)*c + c*im*K
    end
    e1 = rand(dist_energy)
    #el = rand(dist_energy,c)
    K = K_sym(e1, beta, el)
    omega_cm1 = im*(lambda-epsilon2*im)*exp(beta*e1)*c + sum(im*obar.*K./(im*K .+ obar))
    res[1] = real(exp(beta*e1)*c/omega_cm1)
    
    for j in 2:ensemble
        for i in 1:c
            e1 = rand(dist_energy)
            el[i] = e1
            K = meanK(e1, beta)
            obar[i] =  im*(lambda-epsilon*im)*exp(beta*e1)*c + c*im*K
        end
        e1 = rand(dist_energy)
       # el = rand(dist_energy,c)
        K = K_sym(e1, beta, el)
        omega_cm1 = im*(lambda-epsilon2*im)*exp(beta*e1)*c + sum(im*obar.*K./(im*K .+ obar))
        res[j] = real(exp(beta*e1)*c/omega_cm1)
    end
    mean(res)*1/pi
end

function omegabar(lambda,c, ep)
    l = lambda -ep*im
    1/4*(im*(-2+c+2*c*l)+sqrt(-4+4*c-c^2-4*c^2*l - 4*c^2*l^2))
end

function rho_barrat_tinfty(lambda, c, T, ensemble, epsilon2)
    beta = 1.0/T
    res = zeros(ensemble);
    dist_energy = Exponential()
    
    obar = omegabar(lambda, c, epsilon2)
    energies = rand(dist_energy,c)
    e1 = rand(dist_energy)
    Kij = K_sym(e1, beta, energies)
    omega_cm1 =  im*(lambda-epsilon2*im)*exp(beta*e1)*c + sum(im*(Kij.*obar)./(im*(Kij) .+ obar))
    res[1] = real(exp(beta*e1)*c/omega_cm1)
    for j in 2:ensemble
        energies = rand(dist_energy,c)
        e1 = rand(dist_energy)
        Kij = K_sym(e1, beta, energies)
        omega_cm1 =  im*(lambda-epsilon2*im)*exp(beta*e1)*c + sum(im*(Kij.*obar)./(im*(Kij) .+ obar))
        res[j] = real(exp(beta*e1)*c/omega_cm1)
    end
    mean(res)*1/pi
end

function rho_barrat_2nd_tinfty(lambda, c, T, Np, ensemble, nsteps, epsilon2)
    beta = 1.0/T
    poparray = init_array(Np)
    res = zeros(ensemble);
    dist_energy = Exponential()
    obar = omegabar(lambda, c, epsilon2)
    os = zeros(Complex{Float64}, c)
    e1_shell = zeros(c)
    for i in 1:c
        energies = rand(dist_energy, c-1)
        e1 = rand(dist_energy)
        e1_shell[i] = e1
        Kij = K_sym(e1, beta, energies)
        os[i] =  im*(lambda-epsilon2*im)*exp(beta*e1)*c + sum(im*(Kij.*obar)./(im*(Kij) .+ obar))
    end
     #energies = rand(dist_energy, c)
     e1 = rand(dist_energy)
     Kij = K_sym(e1, beta, e1_shell)
    omega_cm1 =  im*(lambda-epsilon2*im)*exp(beta*e1)*c + sum(im*(Kij.*os)./(im*(Kij) .+ os))
    
    res[1] = real(exp(beta*e1)*c/omega_cm1)
    for j in 2:ensemble
        for i in 1:c
            energies = rand(dist_energy, c-1)
            e1 = rand(dist_energy)
            e1_shell[i] = e1
            Kij = K_sym(e1, beta, energies)
            os[i] =  im*(lambda-epsilon2*im)*exp(beta*e1)*c + sum(im*(Kij.*obar)./(im*(Kij) .+ obar))
        end
         #energies = rand(dist_energy, c)
         e1 = rand(dist_energy)
         Kij = K_sym(e1, beta, e1_shell)
        omega_cm1 =  im*(lambda-epsilon2*im)*exp(beta*e1)*c + sum(im*(Kij.*os)./(im*(Kij) .+ os))
        
        res[j] = real(exp(beta*e1)*c/omega_cm1)
    end
    mean(res)*1/pi
end

function omegabarbouchaud(lambda,c, ep)
    l = lambda -ep*im
    1/2*(im*(-2+c+c*l)+sqrt(-4+4*c-c^2-2*c^2*l - c^2*l^2))
end

function rho_bouchaud_tinfty(lambda, c, T, Np, ensemble, nsteps, epsilon)
    tau_dist = Pareto(T)
    res = zeros(ensemble);
    tau = rand(tau_dist)
    obar = omegabarbouchaud(lambda, c, epsilon)
    omega_cm1 = im*(lambda-epsilon*im)*tau*c + c*im*obar./(1.0*im .+ obar)
    res[1] = real(tau*c/omega_cm1)
    for j in 2:ensemble
        tau = rand(tau_dist)
        obar = omegabarbouchaud(lambda, c, epsilon)
        omega_cm1 = im*(lambda-epsilon*im)*tau*c + c*im*obar./(1.0*im .+ obar)
        res[j] = real(tau*c/omega_cm1)
    end
    mean(res)*1/pi
end

function rho_bouchaud_2ndtinfty(lambda, c, T, Np, ensemble, nsteps, epsilon)
    tau_dist = Pareto(T)
    res = zeros(ensemble);
    os = zeros(Complex{Float64}, c)
    obar = omegabar(lambda, c, epsilon)
    for i in 1:c
        tau = rand(tau_dist)
        os[i] = im*(lambda-epsilon*im)*tau*c + (c-1)*im*obar./(1.0*im .+ obar)
    end
    tau = rand(tau_dist)
    omega_cm1 = im*(lambda-epsilon*im)*tau*c + sum(im*os./(1.0*im .+ os))
    res[1] = real(tau*c/omega_cm1)
    
    for j in 2:ensemble
        for i in 1:c
            tau = rand(tau_dist)
            os[i] = im*(lambda-epsilon*im)*tau*c + (c-1)*im*obar./(1.0*im .+ obar)
        end
        tau = rand(tau_dist)
        omega_cm1 = im*(lambda-epsilon*im)*tau*c + sum(im*os./(1.0*im .+ os))
        res[j] = real(tau*c/omega_cm1)
    end
    mean(res)*1/pi
end
