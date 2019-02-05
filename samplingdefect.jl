##Generate a tau variable depending on temperature
function tauj(T)
    dist = Pareto(T)
    rand(dist)
end

##Complete y_j (as arguments it requires tauj, tau and c)
function yj(tj, t, c)
    (1 + 2*t/(c-2) + t/tj )^(-1)
end

##yj under the approximation tau >> c-2 (small lambda). This will yield tilde{Y} value
function yj(tj, c)
    (2/(c-2)+1/tj)^-1
end

##Returns a sample of the complete yj (it requires a tau, temperature and connectivity)
function yjsample(t, T, c )
    tj = tauj(T)
    return yj(tj, t, c)
end

##Returns a sample of the incomplete yj (it requires temperature and c)
function yjsample(T, c)
    tj = tauj(T)
    return yj(tj, c)
end

##It retuns n values of Y, each defined as the sum of the complete y_js over c. They are conditioned on the value of tau
function conditionalY(tau, c, T, n)
    k = 1
    ys = zeros(n, c)
    for i in 1:n
        for j in 1:c
            ys[i,j] = yjsample(tau, T, c)
        end
    end
    Y = [sum(ys[i,1:c]) for i in 1:n]/c;
end

##It retuns n values of Y, each defined as the sum of the incomplete y_js over c. They don't require tau as an argument
function conditionalY(c, T, n)
    k = 1
    ys = zeros(n, c)
    for i in 1:n
        for j in 1:c
            ys[i,j] = yjsample(T, c)
        end
    end
    Y = [sum(ys[i,1:c]) for i in 1:n]/c;
end

##Prefactors based on the correction to the mean field
function prefactors_tlow(l,c,T)
    T*l^(T-1)*(sin(pi*T)/(pi*T))^(T)*((c-2)/2)^(T^2-T)
end

function prefactors_thigh(l,c,T)
    T*l^(T-1)*(c/(c-2))^T
end

