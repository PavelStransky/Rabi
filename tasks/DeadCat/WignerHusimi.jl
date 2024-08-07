using Distributed

@everywhere WORKERS = 8
include("../../Calculation.jl")
@everywhere const PATH = "d:/results/Rabi/deadcat/f2/"
@everywhere pyplot(size = (1200, 1200))

function Run()
    R = 100
    λ = 0.75
    μ = 0.130 / R
    δ = 0.5
    ν = μ
    n = 500
    num = 401

    rabi = Rabi(N=n, R=R, δ=δ, λ=λ, μ=μ, ν=ν)

    # Final values for the figure
    ts = [0, 3, 7, 10, 15, 20]

    # Ground state
    rabi0 = Rabi(N=n, R=R)
    Ψ0 = eigenstates(rabi0)
    Ψ0 = Ψ0[2][1]

    start = 1
    finish = length(ts)

    # wigner = pmap((i)->Wigner(rabi; Ψ0=Ψ0, operators=[:P=>PGS(rabi), :Π=>Parity(rabi), :q=>X(rabi), :Jx=>Jx(rabi)], ts=ts, index=[i], xs=LinRange(-1.7, 1.7, num), ys=LinRange(-0.8, 0.8, num), clim=(-0.1, 0.1), saveData=true, showGraph=false, marginals=false), start:finish)

    rabi = Rabi(R=R, δ=δ, λ=λ, μ=μ, ν=ν)
    wigner = pmap((i)->Wigner(rabi; operators=[:P=>PGS(rabi), :Π=>Parity(rabi), :q=>X(rabi), :Jx=>Jx(rabi)], ts=ts, index=[i], xs=LinRange(-1.7, 1.7, 3*num), ys=LinRange(-0.8, 0.8, 3*num), clim=(-0.1, 0.1), saveData=true, showGraph=false, marginals=true), start:finish)
end

Run()
exit()