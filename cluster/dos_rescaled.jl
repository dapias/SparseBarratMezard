@everywhere include("../src/SparseBarratMezard.jl")
@everywhere using .SparseBarratMezard
@everywhere using JLD


@everywhere function dos(lambda, T; N = 10^4, c = 50,
             ensemble = 10^7, ep2 = [1e-3, 1e-4, 1e-5])
    a = load("../data/populations_rescaled/T=$(T)/lambda=$(lambda)N=$(N)c=$(c)epsi_trescientos.jld")
    pop = a["pop"]
    
    dos = DOSr(lambda, c, T, N, ensemble, ep2, pop)

    return dos
end





@everywhere l_grid = collect(-3:0.02:0.)
#@everywhere l_grid = -exp10.(range(-3,stop =0, length = 101))
@everywhere T = 0.2


dos_pp = pmap(x -> dos(x,T), l_grid)
save("../data/dos_rescaled/T=$(T)c=50.jld", "dos", dos_pp, "lambda",l_grid,
     "ensemble", 10^7, "epsilon2", [1e-3, 1e-4, 1e-5])
#end



    
