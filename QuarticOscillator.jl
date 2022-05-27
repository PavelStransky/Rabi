using QuantumOptics
using Plots

include("Quantum.jl")
include("HOWaveFunction.jl")

struct CUSP <: QuantumSystem
    N
    ħ
    a
end

# If we have the Rabi structure and no new parameters, we copy the values
function Copy(cusp::CUSP; kwargs...)
    N = haskey(kwargs, :N) ? kwargs[:N] : rabi.N
    ħ = haskey(kwargs, :ħ) ? kwargs[:ħ] : rabi.ħ
    a = haskey(kwargs, :a) ? kwargs[:a] : rabi.a

    return CUSP(N, ħ, a)
end

Hamiltonian(p, x; ħ=0.1, a=0.01) = ħ^2 * p^2 / 2 - 2 * x^2 + x^4 + a * x

Basis(cusp::CUSP) = QuantumOptics.FockBasis(cusp.N)

X(cusp::CUSP) = sqrt(cusp.hbar / 2.0) * (destroy(Basis(cusp)) + create(Basis(cusp)))
P(cusp::CUSP) = im * sqrt(cusp.hbar / 2.0) * (create(Basis(cusp)) - destroy(Basis(cusp)))

H(cusp::CUSP) = Hamiltonian(P(cusp), X(cusp); ħ=cusp.ħ, a=cusp.a)

eigenstates(cusp::CUSP) = QuantumOptics.eigenstates(dense(H(cusp)))
eigenstates(cusp::CUSP, i) = QuantumOptics.eigenstates(dense(H(cusp)), i)

function States()
    cusp = CUSP(500, 0.2, 0.0)
    energies, states = eigenstates(cusp)

    xpoints = LinRange(-2, 2, 500)
    p = plot(xpoints, Hamiltonian.(0, xpoints), lw=2)

    for (i, state) in enumerate(states[1:50])
        ypoints = [cusp.hbar^2 * real(WaveFunction(state, x; s=sqrt(1.0 / cusp.hbar))) + energies[i] for x in xpoints]
        p = plot!(p, xpoints, ypoints, ylims=(-1.2, 2))
    end

    display(p)

    return energies
end

function EvolutionLines()
    cusp = CUSP(500, 0.1, 0.01)

    Ψ0 = coherentstate(Basis(cusp), 0)

    mint = 0.0
    maxt = 25.0
    numt = 50

    h = H(cusp)
    u = exp(dense(-im / cusp.hbar * h * mint))
    du = exp(dense(-im / cusp.hbar * h * (maxt - mint) / numt))

    xpoints = LinRange(-2, 2, 1000)
    p = plot(xpoints, Hamiltonian.(0, xpoints), lw=2)

    Ψ = u * Ψ0      # Initial state

    for i = 1:numt
        println(i)
        ypoints = [cusp.hbar * cusp.hbar * (0.5 * real(WaveFunction(Ψ, x; s=sqrt(1.0 / cusp.hbar))) + i - 1) for x in xpoints]
        p = plot!(p, xpoints, ypoints, ylims=(-0.1, 1), legend=false)

        Ψ = du * Ψ
    end

    display(p)
end

function EvolutionHM()
    cusp = CUSP(1000, 0.1, 0.01)

    Ψ0 = coherentstate(Basis(cusp), 0)

    mint = 0.0
    maxt = 50.0
    numt = 500

    h = H(cusp)
    u = exp(dense(-im / cusp.hbar * h * mint))
    du = exp(dense(-im / cusp.hbar * h * (maxt - mint) / numt))

    numx = 601
    xpoints = LinRange(-1.5, 1.5, numx)
    result = Matrix{ComplexF64}(undef, numx, numt)

    Ψ = u * Ψ0      # Initial state

    s = sqrt(1.0 / cusp.hbar)

    for i = 1:numt
        println(i)
        for j = 1:numx
            result[j, i] = WaveFunction(Ψ, xpoints[j]; s=s)
        end
        Ψ = du * Ψ
    end

    p = heatmap(real(result), color=:nipy_spectral, grid=false)
    savefig(p, "real.png")

    p = heatmap(imag(result), color=:nipy_spectral, grid=false)
    savefig(p, "imag.png")
end

pyplot(size=(1200, 1000))

