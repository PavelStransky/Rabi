using QuantumOptics
using Plots

include("Quantum.jl")
include("HOWaveFunction.jl")

pyplot(size=(1200, 1000))

struct SNAIL <: QuantumSystem
    N
    Δ
    K
    ϵ2

    function SNAIL(; N=100, Δ=0, K=1, ϵ2=0)
        return new(N, Δ, K, ϵ2)
    end
end

# If we have the SNAIL structure and no new parameters, we copy the values
function Copy(snail::SNAIL; kwargs...)
    N = haskey(kwargs, :N) && kwargs[:N] > 0 ? kwargs[:N] : snail.N
    Δ = haskey(kwargs, :Δ) ? kwargs[:Δ] : snail.Δ
    K = haskey(kwargs, :K) ? kwargs[:K] : snail.K
    ϵ2 = haskey(kwargs, :ϵ2) ? kwargs[:ϵ2] : snail.ϵ2

    return SNAIL(N=N, Δ=Δ, K=K, ϵ2=ϵ2)
end

Basis(snail::SNAIL) = QuantumOptics.FockBasis(snail.N)

" Operators "
A(snail::SNAIL) = destroy(Basis(snail))
Ad(snail::SNAIL) = create(Basis(snail))
N(snail::SNAIL) = number(Basis(snail))

X(snail::SNAIL) = sqrt(1 / 2.0) * (destroy(Basis(snail)) + create(Basis(snail)))
P(snail::SNAIL) = im * sqrt(1 / 2.0) * (create(Basis(snail)) - destroy(Basis(snail)))

ΨGS(snail::SNAIL) = coherentstate(Basis(snail), 0)

Hamiltonian(p, x, Δ, K, ϵ2) = -Δ / 2 * (x^2 + p^2) + K / 4 * (x^2 + p^2)^2 - ϵ2 * (x^2 - p^2)

AllOperators(snail::SNAIL) = [:q=>X(snail), :p=>P(snail), :n=>N(snail)]

H(snail::SNAIL) = -snail.Δ * create(Basis(snail)) * destroy(Basis(snail)) + snail.K / 4 * create(Basis(snail)) * create(Basis(snail)) * destroy(Basis(snail)) * destroy(Basis(snail)) - snail.ϵ2 * (create(Basis(snail)) * create(Basis(snail)) + destroy(Basis(snail)) * destroy(Basis(snail)))

eigenstates(snail::SNAIL) = QuantumOptics.eigenstates(dense(H(snail)))
eigenstates(snail::SNAIL, i) = QuantumOptics.eigenstates(dense(H(snail)), i)

Size(snail::SNAIL) = 1
Dimension(snail::SNAIL) = snail.N

" Evolution operator "
U(snail::SNAIL, t::Float64) = exp(dense(-im * t * H(snail)))

