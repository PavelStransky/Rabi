using QuantumOptics
using Statistics
using LinearAlgebra

include("Export.jl")
include("Quantum.jl")

pyplot()

" Lipkin type "
struct Lipkin <: QuantumSystem
    j
    ω
    λ

    function Lipkin(; j::Rational=1//2, ω=1, λ=0)
        return new(j, ω, λ)
    end

    function Lipkin(; j::Int=1, ω=1, λ=0)
        return new(j, ω, λ)
    end
end

# If we have the Lipkin structure and no new parameters, we copy the values
function Copy(lipkin::Lipkin; kwargs...)
    j = haskey(kwargs, :j) ? kwargs[:j] : lipkin.j
    ω = haskey(kwargs, :ω) ? kwargs[:ω] : lipkin.ω
    λ = haskey(kwargs, :λ) ? kwargs[:λ] : lipkin.λ

    if haskey(kwargs, :ξ)
        ω = kwargs[:ξ]
        λ = 1 - ω
    end

    return Lipkin(; j=j, ω=ω, λ=λ)
end
    
Size(lipkin::Lipkin) = 2 * lipkin.j
Dimension(lipkin::Lipkin) = Size(lipkin) + 1

""" Print """
Base.String(lipkin::Lipkin) = "Lipkin(j=$(lipkin.j), ω=$(lipkin.ω), λ=$(lipkin.λ))"
Base.show(io::IO, lipkin::Lipkin) = print(io, String(lipkin))

SpinBasis(lipkin::Lipkin) = QuantumOptics.SpinBasis(lipkin.j)

" Operators "
Jx(lipkin::Lipkin) = 0.5 * sigmax(SpinBasis(lipkin))
Jy(lipkin::Lipkin) = 0.5 * sigmay(SpinBasis(lipkin))
Jz(lipkin::Lipkin) = 0.5 * sigmaz(SpinBasis(lipkin))

Jp(lipkin::Lipkin) = 0.5 * sigmap(SpinBasis(lipkin))
Jm(lipkin::Lipkin) = 0.5 * sigmam(SpinBasis(lipkin))

Id(lipkin::Lipkin) = one(SpinBasis(lipkin))

" Hamiltonian "
H0(lipkin::Lipkin) = lipkin.ω * Jz(lipkin)
Hλ(lipkin::Lipkin) = -1 / Size(lipkin) * Jx(lipkin) * Jx(lipkin)

# It is not clear to me why the factor 0.5 is here 
# (but it makes it consistent with all the results from elsewhere)
H(lipkin::Lipkin) = H0(lipkin) + lipkin.λ * Hλ(lipkin)

" Evolution operator "
U(lipkin::Lipkin, t::Float64) = exp(dense(-im * t * H(lipkin)))

" Ground state "
ΨGS(lipkin::Lipkin) = spindown(SpinBasis(lipkin))
DensityMatrix(_::Lipkin, Ψ) = projector(Ψ)

AllOperators(lipkin::Lipkin) = [:Jx=>Jx(lipkin), :Jy=>Jy(lipkin), :Jz=>Jz(lipkin), :P=>PGS(lipkin), :r_r=>PGS(lipkin)]
    
function PlotΨ(lipkin::Lipkin, Ψ)
    result = Matrix{Float64}(undef, 2, Dimension(lipkin))
    m = zeros(Dimension(lipkin))

    if -minimum(real(Ψ.data)) > maximum(real(Ψ.data)) Ψ = -Ψ end

    for i in 1:Dimension(lipkin)
        m[i] = 1
        if i > 1 m[i - 1] = 0 end
        up = Ket(SpinBasis(lipkin), [1, 0])
        down = Ket(SpinBasis(lipkin), [0, 1])

        result[1, i] = real(dagger(down) * Ψ)
        result[2, i] = real(dagger(up) * Ψ)
    end

    p = heatmap(result, color=:bwr, clim=(-0.25, 0.25), yticks=([1.0, 2.0], ["down", "up"]))
    display(p)

    return result, p
end

function eigenstates(lipkin::Lipkin) 
    energies, vectors = QuantumOptics.eigenstates(dense(H(lipkin)))
    return energies / lipkin.j, vectors
end

function eigenstates(lipkin::Lipkin, limit) 
    energies, vectors = QuantumOptics.eigenstates(dense(0.5 * (H(lipkin) + dagger(H(lipkin)))), limit)
    return energies / lipkin.j, vectors
end
