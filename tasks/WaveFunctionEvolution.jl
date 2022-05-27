const PATH = ""

include("../Rabi.jl")
include("../HOWaveFunction.jl")

pyplot(size=(1000, 800))

function WaveFunctionEvolution()
    rabi = Rabi(N=300, R=100, λ=0.75, δ=0.0, μ=0.001, ν=0.001)
    Ψ0 = ΨGS(rabi)

    ts = LinRange(0, 25, 51)

    print("Schrödinger ", rabi, "...")
    time = @elapsed tout, Ψt = timeevolution.schroedinger(ts, Ψ0, H(rabi))
    println(time, "s")

    x = LinRange(-2.0, 2.0, 2001)
    for (i, Ψ) in enumerate(Ψt)
        tstr = format(tout[i], precision=2)
        print("T[$i]=$tstr...")

        ρ = ptrace(Ψ, 2)
        time = @elapsed y = [WaveFunctionSquare(ρ, x1, s=sqrt(2.0 * rabi.R)) for x1 in x]
        println(time)

        p = plot(x, y, legend=false, title="\$t=$tstr\$")
        display(p)
        savefig(p, "$(PATH)WF_$(rabi)_$i.png")
    end
end