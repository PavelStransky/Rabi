include("../Rabi.jl")

const PATH = "d:/results/rabi/hrabac/"

pyplot(size=(800, 600))

""" Survival probabilities P, Pb, Pq """
function SurvivalProbabilities(rabi::Rabi)
    Ψ0 = ΨGS(rabi)
    len = 1001
    ts = LinRange(0, 100, len)

    # Schroedinger time evolution
    print("Schrödinger ", rabi, "...")
    time = @elapsed tout, Ψt = timeevolution.schroedinger(ts, Ψ0, H(rabi))
    println(time, "s")

    pb = Array{Float64}(undef, len)
    pq = Array{Float64}(undef, len)
    p = Array{Float64}(undef, len)

    Ψ0q = spindown(SpinBasis(rabi))
    Ψ0b = coherentstate(FockBasis(rabi), 0)
    Ψ0 = tensor(Ψ0b, Ψ0q)

    for (i, Ψ) in enumerate(Ψt)
        ρ = projector(Ψ)
        ρq = ptrace(Ψ, 1)
        ρb = ptrace(Ψ, 2)

        pb[i] = dagger(Ψ0b) * ρb * Ψ0b
        pq[i] = dagger(Ψ0q) * ρq * Ψ0q

        p[i] = dagger(Ψ0) * ρ * Ψ0
    end

    pl = plot(ts, p)
    pl = plot!(pl, ts, pq)
    pl = plot!(pl, ts, pb)
    display(plot(pl))

    return p, pb, pq
end

#SurvivalProbabilities(Rabi(N=500, R=100, λ=0.75, δ=0.5, μ=0.4, ν=0.4))

function SurvivalProbabilitiesConvergence(rabi::Rabi)
    minR = 1
    maxR = 3
    numR = 21

    Rs = 10 .^ LinRange(minR, maxR, numR)

    jx = Array{Float64}(undef, numR)
    jy = Array{Float64}(undef, numR)
    jz = Array{Float64}(undef, numR)
    p = Array{Float64}(undef, numR)

    for n in 1:numR
        R = (trunc(Int, Rs[n]) ÷ 2) * 2

        rabi = Copy(rabi; R=R, N=-1)
        rabi = Copy(rabi; N=2 * rabi.N)

        Ψ0 = Ket(FockBasis(rabi))
        v = zeros(length(Ψ0.data))
        v[R ÷ 2 + 1] = 1.0
        Ψ0.data = v
        Ψ0 = tensor(Ψ0, spindown(SpinBasis(rabi)))
        
        operators = [:Jx=>Jx(rabi), :Jy=>Jy(rabi), :Jz=>Jz(rabi), :P=>projector(Ψ0)]

        jx[n], jy[n], jz[n], p[n] = AsymptoticValuesMatrix(rabi, operators; Ψ0=Ψ0)
    end    

    return Rs, jx, jy, jz, p
end

Rs, jx, jy, jz, p = SurvivalProbabilitiesConvergence(Rabi(δ=0.5, λ=0.75))
display(plot(Rs, jz, xaxis=:log, xlabel="R", ylabel="Jz", legend="State E=0 (R/2 + 1)"))