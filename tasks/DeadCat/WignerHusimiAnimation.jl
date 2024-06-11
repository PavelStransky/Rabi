using Distributed

@everywhere WORKERS = 8
include("../../Calculation.jl")
@everywhere const PATH = "d:/results/Rabi/deadcat/"
@everywhere pyplot(size = (1200, 1200))

function Run()
    R = 100
    λ = 0.75
    μ = 0.135341379 / R
    δ = 0.0
    ν = μ
    n = 501
    num = 1000

    rabi = Rabi(N=n, R=R, δ=δ, λ=λ, μ=μ, ν=ν)

    ts = LinRange(0, 100, 2001)
    #ts = LinRange(0, 50, 101)
    #ts = LinRange(0, 30, 151)

    #Squeezed state
    # Ψ0 = tensor(squeeze(FockBasis(rabi), 1) * coherentstate(FockBasis(rabi), 0), spindown(SpinBasis(rabi)))
    
    # Ground state
    rabi0 = Rabi(N=n, R=R)
    Ψ0 = eigenstates(rabi0)
    Ψ0 = Ψ0[2][1]

    start = 451
    finish = length(ts)
    finish = 500

    wigner = pmap((i)->Wigner(rabi; Ψ0=Ψ0, operators=[:P=>PGS(rabi), :Π=>Parity(rabi), :q=>X(rabi), :Jx=>Jx(rabi)], ts=ts, index=[i], xs=LinRange(-1.4, 1.4, num), ys=LinRange(-0.7, 0.7, num), clim=(-0.1, 0.1), saveData=false, showGraph=false, marginals=true), start:finish)
    husimi = pmap((i)->Wigner(rabi; Ψ0=Ψ0, husimi=true, operators=[:P=>PGS(rabi), :Π=>Parity(rabi), :q=>X(rabi), :Jx=>Jx(rabi)], ts=ts, index=[i], xs=LinRange(-1.4, 1.4, num), ys=LinRange(-0.7, 0.7, num), clim=(-0.04, 0.04), saveData=false, showGraph=false, marginals=true), start:finish)
end

Run()
exit()