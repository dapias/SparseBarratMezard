using Distributions

function generate_poisson_matrix(N, c)
    values = [1.,-1.,0.]
    probabilities = [1/2*c/N, 1/2*c/N, (1- c/N)]
    dist = Categorical(probabilities)
    
    A = Symmetric(values[rand(dist,N,N)]);
    for j in 1:N
        A[j,j] = 0.
    end
    A
end

function generate_gaussian_matrix(N,c)
    A = zeros(N,N)
    for i in 1:N
        for j in 1:N
            if j > i
                if rand() <= c/N
                    A[i,j] = rand(Normal(0.,sqrt(1/c)))
                end
            end
        end
    end

    A = Symmetric(A)
end


function neighbors_list(M)
    v = Array{Array{Int64,1}}(0)
    N = size(M)[1]
    for i in 1:N
        v_i = Array{Int64}(0)
        for j in 1:N
            if M[i,j] != 0.
                push!(v_i, j)
            end
        end
        v = vcat(v, [v_i])    
    end
    v
end

function local_error(lambda, matrix, neighs, Delta_old, epsilon)
    N = size(matrix)[1]
  
    Delta_new = zeros(Complex{Float64}, N,N);
    Delta_sum = zeros(Complex{Float64}, N,N)
    
    ##Sum involved in Δ^(j)_i, i.e. ∑_{k ∈ ∂i \ j} is stored in Delta_sum
    for i in 1:N
        for j in 1:N
            if j!= i
                for k in neighs[i]
                    if k != j
                        Delta_sum[i,j] += matrix[i,k]^2.*Delta_old[k,i]
                    end
                end
             end
        end
    end
        
    for i in 1:N
        for j in 1:N
            if j!=i
                Delta_new[i,j] = 1./(lambda - epsilon - Delta_sum[i,j])
            end
        end
    end

    
    error = sum(abs.(Delta_new .- Delta_old))
        
    return error,Delta_new
end

function rhocavity(x, c, N; epsilon = 0.005im, tolerance = 0.1, generator = generate_poisson_matrix)
    A = generator(N,c)
    nei = neighbors_list(A)
    error = 10.*tolerance*N*N
    
    Delta = zeros(Complex{Float64}, N,N);   
    
       
    while error > tolerance
        error, Delta = local_error(x, A,nei, Delta, epsilon)
    end
    
    sum_var = 0.
    for i in 1:N
        sum_i = 0.
        for k in nei[i]
            sum_i += A[i,k]^2.*Delta[k,i]
        end
        var_i = 1./(x - epsilon - sum_i)
        sum_var += imag(var_i)
    end
   
    return 1/(pi*N)*sum_var
end        
