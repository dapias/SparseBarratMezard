export generate_sparse_barrat_matrix, generate_barrat_matrix, rho, DOSIPR

###Generates the instances
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

function generate_barrat_matrix(n::Int64,c::Int64, beta::Float64)

    Ms, energias, nei = generate_sparse_barrat_matrix(n,c,beta)

    Ms = Matrix(Symmetric(Ms))
    
    return Ms, energias, nei
end

####Cavity system
function symmetric_f(e1::Float64, beta::Float64, e2::Float64)
    exp(beta*(e1 + e2)/2)/(2*cosh(beta*(e1-e2)/2.)) 
end

function local_error(lambda::Float64,  neighs::Array{Array{Int64,1},1},
        Omega_old::SparseMatrixCSC{Complex{Float64},Int64}, epsilon::Float64,  beta::Float64, energies::Array, c::Int64, fs::SparseMatrixCSC{Complex{Float64},Int64}, N::Int64, Omega_sum::SparseMatrixCSC{Complex{Float64},Int64})

    reg = 10^-10.
    
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
###Two ways of computing the error: Absolute or relative. The absolute demands less memory (used here)
    #relative_real = real.(Omega_sum .- Omega_old)./(abs.(real.(Omega_old))  .+ reg)
    #relative_imag = imag.(Omega_sum .- Omega_old)./(abs.(imag.(Omega_old)) .+ reg)
    #max_real = maximum(abs.(relative_real))
    #max_imag = maximum(abs.(relative_imag))
    #error = maximum([max_real, max_imag])
    error = sum(abs.(Omega_sum .- Omega_old))  

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
    

function DOSIPR(lambda::Float64, c::Int64, n::Int64, beta::Float64, epsilon::Float64, A, energies, nei; tolerance = 0.1)
    
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
        sum_ipr += norm(1.0/(omega_j*(exp(-beta*energies[j])/c)))^2
        sum_var += real(1/(omega_j*(exp(-beta*energies[j])/c)))
    end
    
    dos = 1/(pi*n)*sum_var
    ipr = 1/n*sum_ipr
    
    return dos, epsilon*ipr/(dos*pi)
end





