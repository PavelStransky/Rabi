using QuantumOptics
using Plots

include("Quantum.jl")
include("HOWaveFunction.jl")

pyplot(size=(1200, 1000))

struct CUSP <: QuantumSystem
    N
    ħ
    a
    b

    function CUSP(; N=100, ħ=0.1, a=0.0, b=-2.0)
        return new(N, ħ, a, b)
    end
end

# If we have the CUSP structure and no new parameters, we copy the values
function Copy(cusp::CUSP; kwargs...)
    N = haskey(kwargs, :N) && kwargs[:N] > 0 ? kwargs[:N] : cusp.N
    ħ = haskey(kwargs, :ħ) ? kwargs[:ħ] : cusp.ħ
    a = haskey(kwargs, :a) ? kwargs[:a] : cusp.a
    b = haskey(kwargs, :b) ? kwargs[:b] : cusp.b

    return CUSP(N=N, ħ=ħ, a=a, b=b)
end

Hamiltonian(p, x; ħ=0.1, a=0.01, b=-2) = ħ^2 * p^2 / 2 + b * x^2 + x^4 + a * x

AllOperators(cusp::CUSP) = [:q=>X(cusp), :p=>P(cusp)]

Basis(cusp::CUSP) = QuantumOptics.FockBasis(cusp.N)

X(cusp::CUSP) = sqrt(cusp.ħ / 2.0) * (destroy(Basis(cusp)) + create(Basis(cusp)))
P(cusp::CUSP) = im * sqrt(cusp.ħ / 2.0) * (create(Basis(cusp)) - destroy(Basis(cusp)))

H(cusp::CUSP) = Hamiltonian(P(cusp), X(cusp); ħ=cusp.ħ, a=cusp.a, b=cusp.b)

eigenstates(cusp::CUSP) = QuantumOptics.eigenstates(dense(H(cusp)))
eigenstates(cusp::CUSP, i) = QuantumOptics.eigenstates(dense(H(cusp)), i)

Size(cusp::CUSP) = 1 / cusp.ħ
Dimension(cusp::CUSP) = cusp.N

" Evolution operator "
U(cusp::CUSP, t::Float64) = exp(dense(-im * t * H(cusp)))

" Parity operator "
Parity(cusp::CUSP) = Operator(Basis(cusp), Matrix(Diagonal([isodd(i) ? -1.0 : 1.0 for i in 1:(cusp.N+1)])))

function States()
    cusp = CUSP(500, 0.2, 0.0)
    energies, states = eigenstates(cusp)

    xpoints = LinRange(-2, 2, 500)
    p = plot(xpoints, Hamiltonian.(0, xpoints), lw=2)

    for (i, state) in enumerate(states[1:50])
        ypoints = [cusp.ħ^2 * real(WaveFunction(state, x; s=sqrt(1.0 / cusp.ħ))) + energies[i] for x in xpoints]
        p = plot!(p, xpoints, ypoints, ylims=(-1.2, 2))
    end

    display(p)

    return energies
end

function EvolutionLines()
    cusp = CUSP(N=1000, ħ=0.1, a=0.01)

    Ψ0 = coherentstate(Basis(cusp), 0)

    mint = 0.0
    maxt = 25.0
    numt = 50

    h = H(cusp)
    u = exp(dense(-im / cusp.ħ * h * mint))
    du = exp(dense(-im / cusp.ħ * h * (maxt - mint) / numt))

    xpoints = LinRange(-2, 2, 1000)
    p = plot(xpoints, Hamiltonian.(0, xpoints), lw=2)

    Ψ = u * Ψ0      # Initial state

    for i = 1:numt
        println(i)
        ypoints = [cusp.ħ * cusp.ħ * (0.5 * real(WaveFunction(Ψ, x; s=sqrt(1.0 / cusp.ħ))) + i - 1) for x in xpoints]
        p = plot!(p, xpoints, ypoints, ylims=(-0.1, 1), legend=false)

        Ψ = du * Ψ
    end

    display(p)
end

function EvolutionHM()
    cusp = CUSP(N=1000, ħ=0.1, a=0.01)

    Ψ0 = coherentstate(Basis(cusp), 0)

    mint = 0.0
    maxt = 50.0
    numt = 500

    h = H(cusp)
    u = exp(dense(-im / cusp.ħ * h * mint))
    du = exp(dense(-im / cusp.ħ * h * (maxt - mint) / numt))

    numx = 601
    xpoints = LinRange(-1.5, 1.5, numx)
    result = Matrix{ComplexF64}(undef, numx, numt)

    Ψ = u * Ψ0      # Initial state

    s = sqrt(1.0 / cusp.ħ)

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
