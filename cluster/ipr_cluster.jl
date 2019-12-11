@everywhere using JLD
@everywhere include("cavitymethod2.jl")



@everywhere function ipr(lambda, T; N = 10^4, c = 3,
             ensemble = 10^7)
    a = load("localization/scaling/population/T=$(T)/lambda=$(lambda)N=$(N).jld")
    pop = a["pop"]
    epsi = a["epsilon"]

    dos, ip = DOSIPR(lambda, c, T, N, ensemble, [epsi], pop, epsilon = epsi)

    return ip
end

@everywhere l_grid = collect(-1:0.01:0.)
@everywhere T_grid = 0.3:0.02:0.38
println("starting")
for T in T_grid
    println(T)
    ip = pmap(x -> ipr(x,T), l_grid)
    save("localization/scaling/ipr/T=$(T).jld", "ipr", ip, "lambda",l_grid)
end



    
