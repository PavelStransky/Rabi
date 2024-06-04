using QuantumOptics
using Statistics
using LinearAlgebra

include("Export.jl")
include("Quantum.jl")

pyplot()

" Rabi (Extended Dicke) type "
struct Rabi <: QuantumSystem
    N
    j
    ω
    R
    λ
    δ
    μ
    ν
    η

    function Rabi(; N::Int=-1, j::Rational=1//2, ω=1, R=100, δ=0, λ=0, μ=0, ν=0, η=0)
        # Heuristic optimal value of the size parameter N (for the given λ, μ, R)
        if N < 0 N = max(-N, Nopt(R, λ, μ, ν)) end
    
        return new(N, j, ω, R, λ, δ, μ, ν, η)
    end
end

# If we have the Rabi structure and no new parameters, we copy the values
function Copy(rabi::Rabi; kwargs...)
    N = haskey(kwargs, :N) ? kwargs[:N] : rabi.N
    j = haskey(kwargs, :j) ? kwargs[:j] : rabi.j
    ω = haskey(kwargs, :ω) ? kwargs[:ω] : rabi.ω
    R = haskey(kwargs, :R) ? kwargs[:R] : rabi.R
    λ = haskey(kwargs, :λ) ? kwargs[:λ] : rabi.λ
    δ = haskey(kwargs, :δ) ? kwargs[:δ] : rabi.δ
    μ = haskey(kwargs, :μ) ? kwargs[:μ] : rabi.μ
    ν = haskey(kwargs, :ν) ? kwargs[:ν] : rabi.ν
    η = haskey(kwargs, :η) ? kwargs[:η] : rabi.η

    # Heuristic optimal value of the size parameter N (for the given λ, μ, R)
    if N < 0 N = max(-N, Nopt(R, λ, μ, ν)) end

    return Rabi(; N=N, j=j, ω=ω, R=R, λ=λ, δ=δ, μ=μ, ν=ν, η=η)
end
    
Size(rabi::Rabi) = rabi.R
Dimension(rabi::Rabi) = rabi.N + 1

DensityMatrix(_::Rabi, Ψ) = ptrace(Ψ, 2)

""" Print """
Base.String(rabi::Rabi) = "Rabi(N=$(rabi.N), ω=$(rabi.ω), R=$(rabi.R), δ=$(rabi.δ), λ=$(rabi.λ), " * (rabi.μ == rabi.ν ? "μ=ν=$(rabi.μ)" : "μ=$(rabi.μ), ν=$(rabi.ν)") * ", η=$(rabi.η))"
Base.show(io::IO, rabi::Rabi) = print(io, String(rabi))

" Optimal value of N "
Nopt(R, λ=0, μ=0, ν=0) = trunc(Int, max(R, R + 3 * R * (λ^2 + 4*μ)))
Nopt(rabi::Rabi) = Nopt(rabi.R, rabi.λ, rabi.μ)

FockBasis(rabi::Rabi) = QuantumOptics.FockBasis(rabi.N)
SpinBasis(rabi::Rabi) = QuantumOptics.SpinBasis(rabi.j)

" Operators "
A(rabi::Rabi) = tensor(destroy(FockBasis(rabi)), one(SpinBasis(rabi)))
Ad(rabi::Rabi) = tensor(create(FockBasis(rabi)), one(SpinBasis(rabi)))
N(rabi::Rabi) = tensor(number(FockBasis(rabi)), one(SpinBasis(rabi)))

Jx(rabi::Rabi) = rabi.j * tensor(one(FockBasis(rabi)), sigmax(SpinBasis(rabi)))
Jy(rabi::Rabi) = rabi.j * tensor(one(FockBasis(rabi)), sigmay(SpinBasis(rabi)))
Jz(rabi::Rabi) = rabi.j * tensor(one(FockBasis(rabi)), sigmaz(SpinBasis(rabi)))

Jp(rabi::Rabi) = rabi.j * tensor(one(FockBasis(rabi)), sigmap(SpinBasis(rabi)))
Jm(rabi::Rabi) = rabi.j * tensor(one(FockBasis(rabi)), sigmam(SpinBasis(rabi)))

Id(rabi::Rabi) = tensor(one(FockBasis(rabi)), one(SpinBasis(rabi)))
M(rabi::Rabi) = (N(r) + Jz(r) + rabi.j * Id(rabi))

X(rabi::Rabi) = 1.0 / sqrt(2.0 * rabi.R) * (A(rabi) + Ad(rabi))
P(rabi::Rabi) = im / sqrt(2.0 * rabi.R) * (Ad(rabi) - A(rabi))

" Hamiltonian "
H0(rabi::Rabi) = rabi.ω * N(rabi) + rabi.R * rabi.ω * Jz(rabi)
Hλ(rabi::Rabi) = (1 + rabi.δ) * (Ad(rabi) * Jm(rabi) + A(rabi) * Jp(rabi)) + (1 - rabi.δ) * (Ad(rabi) * Jp(rabi) + A(rabi) * Jm(rabi))
Hμ(rabi::Rabi) = (Ad(rabi) + A(rabi)) * Jz(rabi)
Hν(rabi::Rabi) = rabi.j * (Ad(rabi) + A(rabi))
Hη(rabi::Rabi) = Jx(rabi)

# It is not clear to me why the factor 2 is here 
# (but it makes it consistent with Dicke in CollectiveModels and with the Python code)
H(rabi::Rabi) = H0(rabi) + 2.0 * sqrt(rabi.R) * (rabi.λ * Hλ(rabi) + rabi.μ * Hμ(rabi) + rabi.ν * Hν(rabi) + rabi.η * Hη(rabi))

" Evolution operator "
U(rabi::Rabi, t::Float64) = exp(dense(-im * t * H(rabi)))

" Ground state "
ΨGS(rabi::Rabi) = tensor(coherentstate(FockBasis(rabi), 0), spindown(SpinBasis(rabi)))

AllOperators(rabi::Rabi) = [:Jx=>Jx(rabi), :Jy=>Jy(rabi), :Jz=>Jz(rabi), :P=>PGS(rabi), :q=>X(rabi), :p=>P(rabi), :n=>N(rabi), :r_r=>PGS(rabi)]

" Parity violating ground-state (Puebla2020)"
αsp(rabi::Rabi) = sqrt(rabi.R) * sqrt((2 * rabi.λ)^2 - (2 * rabi.λ)^-2) / 2
ssp(rabi::Rabi) = -0.25 * log(1.0 - (2 * rabi.λ)^-4)
Displacement(rabi::Rabi) = exp(dense(αsp(rabi) * create(FockBasis(rabi)) - conj(αsp(rabi)) * destroy(FockBasis(rabi))))
Squeeze(rabi::Rabi) = exp(dense(0.5 * (conj(ssp(rabi)) * create(FockBasis(rabi)) * create(FockBasis(rabi)) - ssp(rabi) * destroy(FockBasis(rabi)) * destroy(FockBasis(rabi)))))
ΨGSP(rabi::Rabi) = tensor(Displacement(rabi) * Squeeze(rabi) * coherentstate(FockBasis(rabi), 0), -sqrt(0.5 * (1 - (2 * rabi.λ)^-2)) * spinup(SpinBasis(rabi)) + sqrt(0.5 * (1 + (2 * rabi.λ)^-2)) * spindown(SpinBasis(rabi)))
PGSP(rabi::Rabi) = projector(ΨGSP(rabi))

" Parity operator "
Parity(rabi::Rabi) = tensor(Operator(FockBasis(rabi), Matrix(Diagonal([isodd(i) ? -1.0 : 1.0 for i in 1:(rabi.N+1)]))), sigmaz(SpinBasis(rabi)))
    
function ProjectParity(rabi::Rabi, vector1::Ket, vector2::Ket)
    p = Parity(rabi)
    m = real([dagger(vector1) * p * vector1 dagger(vector1) * p * vector2;
        dagger(vector2) * p * vector1 dagger(vector2) * p * vector2])
    
    println(m)
    ev = eigen(m)
    println(ev)

    a = vector1 * ev.vectors[1, 1] + vector2 * ev.vectors[2, 1]
    b = vector1 * ev.vectors[1, 2] + vector2 * ev.vectors[2, 2]

    return a, b
end

function PlotΨ(rabi::Rabi, Ψ)
    result = Matrix{Float64}(undef, 2, rabi.N + 1)
    m = zeros(rabi.N + 1)

    if -minimum(real(Ψ.data)) > maximum(real(Ψ.data)) Ψ = -Ψ end

    for i in 1:(rabi.N + 1)
        m[i] = 1
        if i > 1 m[i - 1] = 0 end
        up = tensor(Ket(FockBasis(rabi), m), Ket(SpinBasis(rabi), [1, 0]))
        down = tensor(Ket(FockBasis(rabi), m), Ket(SpinBasis(rabi), [0, 1]))

        result[1, i] = real(dagger(down) * Ψ)
        result[2, i] = real(dagger(up) * Ψ)
    end

    p = heatmap(result, color=:bwr, clim=(-0.25, 0.25), yticks=([1.0, 2.0], ["down", "up"]))
    display(p)

    return result, p
end

function eigenstates(rabi::Rabi) 
    energies, vectors = QuantumOptics.eigenstates(dense(H(rabi)))
    return energies / rabi.R, vectors
end

function eigenstates(rabi::Rabi, limit) 
    energies, vectors = QuantumOptics.eigenstates(dense(H(rabi)), limit)
    return energies / rabi.R, vectors
end
