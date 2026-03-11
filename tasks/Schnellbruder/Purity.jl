include("../../Rabi.jl")

const PATH = "d:/results/rabi/Schnellbruder/"

function CalculatePurity(rabii, λf; maxt=200, numt=2001)
    # Initial and final systems
    rabif = Copy(rabii; λ=λf)
    println(rabii.N)

    Ψ0 = SingleWellState(rabii)	

    ts = LinRange(0.0, maxt, numt)
    purity = Array{Any}(undef, numt)

    time = @elapsed tout, Ψt = timeevolution.schroedinger(ts, Ψ0, H(rabif))
    println(time, "s")

    for (i, Ψ) in enumerate(Ψt)
        ρq = ptrace(Ψ, 2)
        purity[i] = real(tr(ρq * ρq))
    end

    # ev = ExpectationValues(rabif, [:Jx=>Jx(rabif)*Jx(rabif), :Jy=>Jy(rabif)*Jy(rabif), :Jz=>Jz(rabif)*Jz(rabif)]; Ψ0=Ψ0, maxt=maxt, numt=2001, asymptotics=false, showGraph=false, saveGraph=false)
    ev = ExpectationValues(rabif, [:Jx=>Jx(rabif), :Jy=>Jy(rabif), :Jz=>Jz(rabif)]; Ψ0=Ψ0, maxt=maxt, numt=2001, asymptotics=false, showGraph=false, saveGraph=false)
    jx = ev[1, :]
    jy = ev[2, :]
    jz = ev[3, :]
    
    p = plot(ts, purity, xlabel="\$t\$", ylabel="Purity", title="Evolution of the purity", legend=false)
    savefig(p, "$(PATH)purity_$(rabii)_$(λf).png")

    Export("$(PATH)purity_ptrace_$(rabii)_$(λf)", tout, purity)

    x = 1 / (2 * rabif.j + 1)
    Export("$(PATH)purity_j_$(rabii)_$(λf)", tout, (jx .* jx .+ jy .* jy .+ jz .* jz) ./ (rabif.j)^2 .* (1 - x) .+ x)
end

CalculatePurity(Rabi(R=50, λ=1.5, δ=0.0, j=1//2), 0.2, maxt=300)