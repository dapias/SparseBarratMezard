using Distributions, LightGraphs, LinearAlgebra

function generate_bouchaud_matrix(n,c, beta)
    dene = Exponential()
    trans_rate(b, con, energ, grafo) = [exp(-b*energ[j])/con for i in vertices(grafo), j in vertices(grafo)]
    L = random_regular_graph(n,c)
    adj = adjacency_matrix(L);
    energias = rand(dene,nv(L))
    rates = trans_rate(beta, c, energias, L);
    M = rates.*adj
    for i in vertices(L)
        M[i,i] = -sum(M[:,i])
    end
    p_eq =  exp.(beta*energias)
    p_eq /= sum(p_eq)
    p_inv = p_eq.^(-1/2)
    p_dir = p_eq.^(1/2)
    Ms = Diagonal(p_inv)*M*Diagonal(p_dir)
    Ms = Matrix(Ms)
    return Ms, rates[1,:]
end

transition_rates(b, con, energ, N::Int64) = 1.0/con*[1.0/(1.0+exp(-b*(energ[i] - energ[j]))) for i in 1:N, j in 1:N]

function generate_barrat_matrix(n,c, beta)
    transition_rates(b, con, energ, grafo) = 1.0/con*[1.0/(1.0+exp(-b*(energ[i] - energ[j]))) for i in vertices(grafo), j in vertices(grafo)]
    
    dene = Exponential()
    L = random_regular_graph(n,c)
    adj = adjacency_matrix(L);
    energias = rand(dene,nv(L))
    rates = transition_rates(beta, c, energias, L);
    M = rates.*adj
    for i in vertices(L)
        M[i,i] = -sum(M[:,i])
    end
    peq = exp.(beta*energias)
    peq /= sum(peq)
    p_inv = peq.^(-1/2)
    p_dir = peq.^(1/2)
    Ms = Diagonal(p_inv)*M *Diagonal(p_dir)
    Ms = Matrix(Ms)
    return Ms, energias
end

function neighbors_list(M)
    v = [Int64[]]
    N = size(M)[1]
    for i in 1:N
        v_i = Int64[]
        for j in 1:N
            if M[i,j] != 0.
                push!(v_i, j)
            end
        end
        v = vcat(v, [v_i])    
    end
    v[2:end]
end

function local_error(lambda, matrix, neighs, Omega_old, epsilon, rates)
    ##Éste es para el problema de Bouchaud
    N = size(matrix)[1]
  
    Omega_new = zeros(Complex{Float64}, N,N);
    Omega_sum = zeros(Complex{Float64}, N,N)
    
    ##Sum involved in ω^(j)_i, i.e. ∑_{k ∈ ∂i \ j} is stored in Omega_sum
    for i in 1:N
        for j in 1:N
            if j!= i
                for k in neighs[i]
                    if k != j
                        Omega_sum[i,j] += Omega_old[k,i]/(1.0+Omega_old[k,i])
                    end
                end
             end
        end
    end
        
    for i in 1:N
        for j in 1:N
            if j!=i
                Omega_new[i,j] = (lambda - im*epsilon)/rates[i] + Omega_sum[i,j]
            end
        end
    end

    
    error = sum(abs.(Omega_new .- Omega_old))
        
    return error, Omega_new
end


function local_error(lambda, matrix, neighs, Omega_old, epsilon, f, rates) ##Cuando pasas "f" usa Barrat calculations
    N = size(matrix)[1]
  
    Omega_new = zeros(Complex{Float64}, N,N);
    Omega_sum = zeros(Complex{Float64}, N,N)

    ##Sum involved in Δ^(j)_i, i.e. ∑_{k ∈ ∂i \ j} is stored in Delta_sum
    for i in 1:N
        for j in 1:N
            if j!= i
                for k in neighs[i]
                    if k != j
                        Omega_sum[i,j] += f[i,k]*Omega_old[k,i]/(f[i,k]+Omega_old[k,i])
                    end
                end
             end
        end
    end
        
    for i in 1:N
        for j in 1:N
            if j!=i
                Omega_new[i,j] = (lambda - epsilon)/rates[i] + Omega_sum[i,j]
            end
        end
    end

    
    error = sum(abs.(Omega_new .- Omega_old))
        
    return error, Omega_new
end



function rhocavity_bouchaud(x, c, N, beta; epsilon = 0.005, tolerance = 0.1)
    A, rates = generate_bouchaud_matrix(N,c, beta)
    nei = neighbors_list(A)
    error = 10.0*tolerance*N*N
    
    Omega = zeros(Complex{Float64}, N,N);   
    
    while error > tolerance
        error, Omega = local_error(x, A,nei, Omega, epsilon, rates)
    end
    
    sum_var = 0.
    for i in 1:N
        sum_i = 0.
        for k in nei[i]
            sum_i += Omega[k,i]/(1.0+ Omega[k,i])
        end
        omega_i = (x - im*epsilon)/rates[i] + sum_i
        sum_var += imag(1.0/(omega_i*rates[i]))
    end
   
    return 1/(pi*N)*sum_var
end

function rhocavity_barrat(x, c, N, beta; epsilon = 0.005im, tolerance = 0.1, generator = generate_barrat_matrix)
    A, energies = generator(N,c, beta)
    nei = neighbors_list(A)
    error2 = 10.0*tolerance*N*N
    
    Omega = zeros(Complex{Float64}, N,N);   

    f = [exp(beta*(energies[i] + energies[j])/2)/(2*cosh(beta*(energies[i]-energies[j])/2.)) for i in 1:N, j in 1:N]
    rates = [exp(-beta*energies[i])/c for i in 1:N]
    
    while error2 > tolerance
        error2, Omega = local_error(x, A,nei, Omega, epsilon, f, rates)
    end
 
    sum_var = 0.
    for i in 1:N
        sum_i = 0.
        for k in nei[i]
            sum_i += Omega[k,i]*f[i,k]/( f[i,k] + Omega[k,i])
        end
        omega_i = (x - epsilon)/rates[i] + sum_i
        sum_var += imag(1.0/(omega_i*rates[i]))
    end
   
    return 1/(pi*N)*sum_var
end        


struct Population
    omegas::Array{Complex{Float64},1}
    energies::Array{Float64,1}
end

function K_sym(e_1, beta, es)
    exp.(beta*(e_1 .+ es)/2)./(2*cosh.(beta*(e_1 .- es)/2.))
end


function init_array(Np)
    dist_energy = Exponential()
    Population(rand(Complex{Float64},Np), rand(dist_energy, Np))
end


function population_update(lambda,c,T,Np, nsteps, epsilon)
    beta = 1.0/T
    poparray = init_array(Np)
    dist_energy = Exponential()
    for i in 1:nsteps
        random_elements = rand(1:Np,c-1)
        omegas = poparray.omegas[random_elements]
        energies = poparray.energies[random_elements]
        e1 = rand(dist_energy)
        Kij = K_sym(e1, beta, energies)
        omega_cm1 =  im*(lambda-epsilon*im)*exp(beta*e1)*c + sum(im*(Kij.*omegas)./(im*(Kij) .+ omegas))
        relement = rand(1:Np)
        poparray.omegas[relement] = omega_cm1
        poparray.energies[relement] = e1
    end
    poparray
end


function population_update!(lambda,c,T,Np, nsteps, epsilon, poparray)
    beta = 1.0/T
    dist_energy = Exponential()
    for i in 1:nsteps
        random_elements = rand(1:Np,c-1)
        omegas = poparray.omegas[random_elements]
        energies = poparray.energies[random_elements]
        e1 = rand(dist_energy)
        Kij = K_sym(e1, beta, energies)
        omega_cm1 =  im*(lambda-epsilon*im)*exp(beta*e1)*c + sum(im*(Kij.*omegas)./(im*(Kij) .+ omegas))
        relement = rand(1:Np)
        poparray.omegas[relement] = omega_cm1
        poparray.energies[relement] = e1
    end
    poparray
end

function rho_barrat_population(lambda, c, T, Np, ensemble, nsteps, epsilon, epsilon2::Float64)
    beta = 1.0/T
    poparray = init_array(Np)
    population_update!(lambda, c, T, Np, nsteps, epsilon, poparray);  ##Equilibrium is reached
    res = zeros(ensemble);
    
    dist_energy = Exponential()
    random_elements = rand(1:Np,c)
    omega_sample = poparray.omegas[random_elements]
    energies = poparray.energies[random_elements]  
    e1 = rand(dist_energy)
    Kij = K_sym(e1, beta, energies)
    omega_c =  im*(lambda-epsilon2*im)*exp(beta*e1)*c + sum(im*(Kij.*omega_sample)./(im*(Kij) .+ omega_sample))
    res[1] = real(exp(beta*e1)*c/omega_c)
 
    for j in 2:ensemble
        random_elements = rand(1:Np,c)
        omega_sample = poparray.omegas[random_elements]
        energies = poparray.energies[random_elements]  
        e1 = rand(dist_energy)
        Kij = K_sym(e1, beta, energies)
        omega_c =  im*(lambda-epsilon2*im)*exp(beta*e1)*c + sum(im*(Kij.*omega_sample)./(im*(Kij) .+ omega_sample))   
        res[j] = real(exp(beta*e1)*c/omega_c)
    end
    mean(res)*1/pi
end



function rho_barrat_population(lambda, c, T, Np, ensemble, epsilon, epsilon2, poparray::Population)
    ##For this function a population array is passed that is assumed is already equilibrated
    beta = 1.0/T
    res = zeros(ensemble);  
    dist_energy = Exponential()
    random_elements = rand(1:Np,c)
    omega_sample = poparray.omegas[random_elements]
    energies = poparray.energies[random_elements]  
    e1 = rand(dist_energy)
    Kij = K_sym(e1, beta, energies)
    omega_c =  im*(lambda-epsilon2*im)*exp(beta*e1)*c + sum(im*(Kij.*omega_sample)./(im*(Kij) .+ omega_sample))
    res[1] = real(exp(beta*e1)*c/omega_c)
    ##update the population with new values
    omega_sample = omega_sample[1:end-1]
    energies = energies[1:end-1]
    Kij = K_sym(e1, beta, energies)
    omega_cm1 = im*(lambda-epsilon*im)*exp(beta*e1)*c + sum(im*(Kij.*omega_sample)./(im*(Kij) .+ omega_sample))
    relement = rand(1:Np)
    poparray.omegas[relement] = omega_cm1
    poparray.energies[relement] = e1
    ###################
    
    for j in 2:ensemble
        random_elements = rand(1:Np,c)
        omega_sample = poparray.omegas[random_elements]
        energies = poparray.energies[random_elements]  
        e1 = rand(dist_energy)
        Kij = K_sym(e1, beta, energies)
        omega_c =  im*(lambda-epsilon2*im)*exp(beta*e1)*c + sum(im*(Kij.*omega_sample)./(im*(Kij) .+ omega_sample))   
        res[j] = real(exp(beta*e1)*c/omega_c)
        ##update the population with new values
        omega_sample = omega_sample[1:end-1]
        energies = energies[1:end-1]
        Kij = K_sym(e1, beta, energies)
        omega_cm1 = im*(lambda-epsilon*im)*exp(beta*e1)*c + sum(im*(Kij.*omega_sample)./(im*(Kij) .+ omega_sample))
        relement = rand(1:Np)
        poparray.omegas[relement] = omega_cm1
        poparray.energies[relement] = e1
        ##########################
    end
    mean(res)*1/pi
end


function population_bouchaud(lambda, c, T, Np, nsteps, epsilon)
    poparray = rand(Complex{Float64},Np)
    tau_dist = Pareto(T)
    for i in 1:nsteps
        random_elements = rand(1:Np,c-1)
        omegas = poparray[random_elements]
        tau = rand(tau_dist)
        omega_cm1 = im*(lambda-epsilon*im)*tau*c + sum(im*omegas./(1.0*im .+ omegas))
        poparray[rand(1:Np)] = omega_cm1
    end
    poparray
end

function rho_bouchaud_population(lambda, c, T, Np, ensemble, nsteps, epsilon, epsilon2::Float64)
    omegas = population_bouchaud(lambda, c, T, Np,  nsteps, epsilon);
    tau_dist = Pareto(T)
    res = zeros(ensemble);
    tau = rand(tau_dist)
    omega_sample = omegas[rand(1:Np,c)]
    omega_cm1 = im*(lambda-epsilon2*im)*tau*c + sum(im*omega_sample./(1.0*im .+ omega_sample))
    res[1] = real(tau*c/omega_cm1)
    for j in 2:ensemble
        tau = rand(tau_dist)
        omega_sample = omegas[rand(1:Np,c)]
        omega_cm1 = im*(lambda-epsilon2*im)*tau*c + sum(im*omega_sample./(1.0*im .+ omega_sample))
        res[j] = real(tau*c/omega_cm1)
    end
    mean(res)*1/pi
end

function rho_bouchaud_population(lambda, c, T, Np, ensemble, epsilon, epsilon2, omegas::Vector{Complex{Float64}})
    tau_dist = Pareto(T)
    res = zeros(ensemble);
    tau = rand(tau_dist)
    omega_sample = omegas[rand(1:Np,c)]
    omega_cm1 = im*(lambda-epsilon2*im)*tau*c + sum(im*omega_sample./(1.0*im .+ omega_sample))
    res[1] = real(tau*c/omega_cm1)
    for j in 2:ensemble
        tau = rand(tau_dist)
        omega_sample = omegas[rand(1:Np,c)]
        omega_cm1 = im*(lambda-epsilon2*im)*tau*c + sum(im*omega_sample./(1.0*im .+ omega_sample))
        res[j] = real(tau*c/omega_cm1)
    end
    mean(res)*1/pi
end

