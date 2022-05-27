using Base: vect
using Distributed
using Quadmath

@everywhere const WORKERS = 6

include("Calculation.jl")

pyplot(size=(1000, 1000))

@everywhere const PATH = "d:\\results\\rabi\\tmp\\"
#@everywhere PATH = ""

#BLAS.set_num_threads(16)

rabi = Rabi(N=100, R=100, λ=0.75, δ=0.5)
energies, vectors = eigenstates(rabi)

Ψ = vectors[1]

Ψ0 = ΨGS(rabi)

ts = LinRange(0, 20, 21)
time = @elapsed tout, Ψt = timeevolution.schroedinger(ts, Ψ0, H(rabi))

xbasis = PositionBasis(-sqrt(2.0 * rabi.R), sqrt(2.0 * rabi.R), 200)
pf = transform(xbasis, FockBasis(rabi))
x = samplepoints(xbasis)

for (j, t) in enumerate(ts)
    Ψf = ptrace(Ψt[j], 2)
    Ψx = pf * Ψf * dagger(pf)

    p = plot(title="t = $t")
    for k = 1:100
        p = plot!(p, x, abs2.(Ψx.data[k,:]), legend=false, ylims=(0, 0.02))
    end
    display(p)
end

#y = abs2.(Ψf.data[j,:])