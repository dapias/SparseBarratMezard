export laplaceinverse

function local_error(s::Union{Float64, Complex{Float64}}, neighs::Array{Array{Int64,1},1}, Omega_old::SparseMatrixCSC,  beta::Float64, energies::Array{Float64,1}, c::Int64,  fs::SparseMatrixCSC{Complex{Float64},Int64}, N::Int64, Omega_sum::SparseMatrixCSC)
   
    for k in 1:N
        value = neighs[k]
        for j in value
            for l in value
                if l != j
                    Omega_sum[k, j] +=  fs[k,l]*Omega_old[l, k]/(fs[k,l] +Omega_old[l,k])
                end
            end
            Omega_sum[k,j] += (s)/(exp(-beta*energies[k])/c) 
        end
    end
    
    error = sum(abs.(Omega_sum .- Omega_old))
    @inbounds fill!(Omega_old, zero(Complex{Float64}))
    
    return error, Omega_sum, Omega_old
end


function rho2(s::Union{Float64, Complex{Float64}}, c::Int64, beta::Float64,  A::Union{Matrix, SparseMatrixCSC}, energies::
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
        error2, Omegas, Omegap = local_error(s, nei,  Omegas,  beta, energies, c, fs, n, Omegap)
    end
 
   sum_var = 0.
    for j in 1:n
        sum_j = 0.
        for k in nei[j]
            sum_j += Omegas[k,j]*fs[k,j]/( fs[k,j] + Omegas[k,j])
        end
        omega_j = s + (exp(-beta*energies[j])/c)*sum_j
        sum_var += (1.0/(omega_j))
    end
   
    return 1/(n)*sum_var
end     

function laplaceinverse(n::Int64, c::Int64, T::Float64)
    beta = 1/T
    M, es, nei = generate_sparse_barrat_matrix(n,c, beta);
    Talbot(s -> rho2(s, c, beta, M, es, nei))   ###Return probability to node 1
end
