using Distributed

@everywhere WORKERS = 4
include("../Calculation.jl")
@everywhere const PATH = "d:/results/Rabi/wigner/"

function Run()
    r = Array{Any}(undef, WORKERS)

    R = 100
    λ = 0.75
    μ = 0.001
    ν = μ
    n = 500
    num = 401

    rabi = Rabi(N=n, R=R, λ=λ, μ=μ, ν=ν)

    ts = LinRange(0, 100, 2001)
    #ts = LinRange(0, 50, 101)

    w2 = WORKERS ÷ 2

    Ψ0 = tensor(squeeze(FockBasis(rabi), 1) * coherentstate(FockBasis(rabi), 0), spindown(SpinBasis(rabi)))
    
    #rabi0 = Rabi(N=n, R=R)
    #Ψ0 = eigenstates(rabi0)
    #Ψ0 = Ψ0[2][5]

    for j in 1:w2
        r[2*j - 1] = @spawn Wigner(rabi; Ψ0=Ψ0, operators=[:P=>PGS(rabi), :Π=>Parity(rabi), :q=>X(rabi), :Jx=>Jx(rabi)], ts=ts, index=j:w2:length(ts), xs=LinRange(-1.4, 1.4, num), ys=LinRange(-0.7, 0.7, num), clim=(-0.1, 0.1), saveData=false, showGraph=false)
        r[2*j] = @spawn Wigner(rabi; Ψ0=Ψ0, husimi=true, operators=[:P=>PGS(rabi), :Π=>Parity(rabi), :q=>X(rabi), :Jx=>Jx(rabi)], ts=ts, index=j:w2:length(ts), xs=LinRange(-1.4, 1.4, num), ys=LinRange(-0.7, 0.7, num), clim=(-0.04, 0.04), saveData=false, showGraph=false)
    end

    for j in 1:WORKERS
        fetch(r[j])
    end
end

Run()