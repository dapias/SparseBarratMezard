include("cavitymethod.jl")
using PyCall

@pyimport scipy.special as sp

function meanK(ene, beta)
    exp(ene*beta)*sp.hyp2f1(1.,1/beta, 1+1/beta, -exp(ene*beta))
end

function rho_barrat_cinfty(lambda, c, T, ensemble, epsilon, epsilon2)
    beta = 1./T
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
    beta = 1./T
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
