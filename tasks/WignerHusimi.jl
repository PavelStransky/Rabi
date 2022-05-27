using Distributed

@everywhere WORKERS = 20
include("../Calculation.jl")
@everywhere const PATH = ""

function Run()
    r = Array{Any}(undef, WORKERS)

    R = 100
    λ = 0.75
    μ = 0.001
    ν = μ
    n = 500
    num = 101

    rabi = Rabi(N=n, R=R, λ=λ, μ=μ, ν=ν)

    ts = LinRange(0, 100, 2001)

    w2 = WORKERS ÷ 2

    for j in 1:w2
        r[2*j - 1] = @spawn Wigner(rabi; operators=[:P=>PGS(rabi), :Π=>Parity(rabi), :q=>X(rabi), :Jx=>Jx(rabi)], ts=ts, index=j:w2:length(ts), xs=LinRange(-1.4, 1.4, num), ys=LinRange(-0.7, 0.7, num), clim=(-0.1, 0.1), saveData=false, showGraph=false)
        r[2*j] = @spawn Wigner(rabi; husimi=true, operators=[:P=>PGS(rabi), :Π=>Parity(rabi), :q=>X(rabi), :Jx=>Jx(rabi)], ts=ts, index=j:w2:length(ts), xs=LinRange(-1.4, 1.4, num), ys=LinRange(-0.7, 0.7, num), clim=(-0.04, 0.04), saveData=false, showGraph=false)
    end

    for j in 1:WORKERS
        fetch(r[j])
    end
end

Run()