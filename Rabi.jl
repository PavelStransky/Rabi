using QuantumOptics
using Statistics
using LinearAlgebra


include("Export.jl")
include("Quantum.jl")

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
    
Size(rabi::Rabi) = 2 * rabi.R * rabi.j
Dimension(rabi::Rabi) = rabi.N + 1

DensityMatrix(_::Rabi, Ψ) = ptrace(Ψ, 2)

""" Print """
function Base.String(rabi::Rabi) 
    result = "Rabi(2j=$(Int(2*rabi.j)), ω=$(rabi.ω), R=$(rabi.R), λ=$(rabi.λ)"
    if rabi.δ != 0
        result *= ", δ=$(rabi.δ)"
    end
    if rabi.μ != 0
        if rabi.μ == rabi.ν
            result *= ", μ=ν=$(rabi.μ)"
        else
            result *= ", μ=$(rabi.μ)"
        end
    end
    if rabi.ν != 0 && rabi.μ != rabi.ν
        result *= ", ν=$(rabi.ν)"
    end
    if rabi.η != 0
        result *= ", η=$(rabi.η)"
    end

    return result * ")"
end

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

Jx(rabi::Rabi) = 0.5 * tensor(one(FockBasis(rabi)), sigmax(SpinBasis(rabi)))
Jy(rabi::Rabi) = 0.5 * tensor(one(FockBasis(rabi)), sigmay(SpinBasis(rabi)))
Jz(rabi::Rabi) = 0.5 * tensor(one(FockBasis(rabi)), sigmaz(SpinBasis(rabi)))

Jp(rabi::Rabi) = 0.5 * tensor(one(FockBasis(rabi)), sigmap(SpinBasis(rabi)))
Jm(rabi::Rabi) = 0.5 * tensor(one(FockBasis(rabi)), sigmam(SpinBasis(rabi)))

Id(rabi::Rabi) = tensor(one(FockBasis(rabi)), one(SpinBasis(rabi)))
M(rabi::Rabi) = (N(r) + Jz(r) + rabi.j * Id(rabi))

X(rabi::Rabi) = 1.0 / sqrt(4.0 * rabi.R * rabi.j) * (A(rabi) + Ad(rabi))
P(rabi::Rabi) = im / sqrt(4.0 * rabi.R * rabi.j) * (Ad(rabi) - A(rabi))

" Hamiltonian "
H0(rabi::Rabi) = rabi.ω * N(rabi) + rabi.R * rabi.ω * (Jz(rabi) + rabi.j * Id(rabi))
Hλ(rabi::Rabi) = (1 + rabi.δ) * (Ad(rabi) * Jm(rabi) + A(rabi) * Jp(rabi)) + (1 - rabi.δ) * (Ad(rabi) * Jp(rabi) + A(rabi) * Jm(rabi))
Hμ(rabi::Rabi) = (Ad(rabi) + A(rabi)) * Jz(rabi)
Hν(rabi::Rabi) = rabi.j * (Ad(rabi) + A(rabi))
Hη(rabi::Rabi) = Jx(rabi)

# It is not clear to me why the factor 2 is here 
# (but it makes it consistent with Dicke in CollectiveModels and with the Python code)
# Coupling normalized to give the critical coupling always at λ = 1
H(rabi::Rabi) = H0(rabi) + 2.0 * sqrt(2.0 * rabi.R * rabi.j) * (0.25 / rabi.j * rabi.λ * rabi.ω * Hλ(rabi) + rabi.μ * Hμ(rabi) + rabi.ν * Hν(rabi) + rabi.η * Hη(rabi))

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
Parity(rabi::Rabi) = tensor(Operator(FockBasis(rabi), Matrix(Diagonal([isodd(i) ? -1.0 : 1.0 for i=0:rabi.N]))),
                        Operator(SpinBasis(rabi), Matrix(Diagonal([isodd(i) ? -1.0 : 1.0 for i=(2*rabi.j:-1:0)]))))
#Parity(rabi::Rabi) = tensor(Operator(FockBasis(rabi), Matrix(Diagonal([isodd(i) ? -1.0 : 1.0 for i in 1:(rabi.N+1)]))), sigmaz(SpinBasis(rabi)))
    
function ProjectParity(rabi::Rabi, vector1::Ket, vector2::Ket)
    p = Parity(rabi)
    m = real([dagger(vector1) * p * vector1 dagger(vector1) * p * vector2;
        dagger(vector2) * p * vector1 dagger(vector2) * p * vector2])
    
    ev = eigen(m)

    a = vector1 * ev.vectors[1, 2] + vector2 * ev.vectors[2, 2]
    b = vector1 * ev.vectors[1, 1] + vector2 * ev.vectors[2, 1]

    return a, b

end

" A parity-violating state for a backward quench "
function SingleWellState(rabi::Rabi; left=true)
    _, vs = eigenstates(rabi, 2) 
    a, b = ProjectParity(rabi, vs[1], vs[2])

    gs1 = (a + b) / sqrt(2)
    gs2 = (a - b) / sqrt(2)

    q1 = ExpectationValue("E", X(rabi), gs1, rabi)

    if q1 < 0 && left
        return gs1
    end

    return gs2
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
    return energies / (2 * rabi.j * rabi.R * rabi.ω), vectors
end

function eigenstates(rabi::Rabi, limit) 
    energies, vectors = QuantumOptics.eigenstates(dense(H(rabi)), limit)
    return energies / (2 * rabi.j * rabi.R * rabi.ω), vectors
end

""" Strength function for the Rabi model and its local extrema """
function StrengthFunction(rabii::Rabi; λf=0.5, threshold=1E-15, envelope_window=3, showgraph=true, savedata=false)
    """ Finds all local extrema in the array a """
    function LocalExtremaIndices(a)
        n = length(a)
        extrema_indices = Int[]

        # Include the first...
        push!(extrema_indices, 1)

        for i in 2:(n-1)
            if (a[i] > a[i-1] && a[i] > a[i+1]) || (a[i] < a[i-1] && a[i] < a[i+1])
                push!(extrema_indices, i)
            end
        end
    
        # ...and the last elements as extrema
        push!(extrema_indices, n)

        return extrema_indices
    end    

    # Initial and final systems
    rabif = Copy(rabii; λ=λf)
    gs = SingleWellState(rabii)

    # Strength function
    sf = StrengthFunction(rabif, gs)

    energies = sf[1]
    probabilities = sf[2]

    indices = findall(p -> p > threshold, probabilities)

    energies = energies[indices]
    probabilities = probabilities[indices]

    envelopex = copy(energies)
    envelopey = copy(probabilities)

    for i in 1:length(probabilities)
        k = argmax(probabilities[max(i - envelope_window, 1):i]) + max(i - envelope_window, 1) - 1
        envelopex[i] = energies[k]
        envelopey[i] = probabilities[k]
    end

        # Find indices of the first occurrences of unique y-values
    uniquey = unique(envelopey)

    unique_indices = Int[]
    for y in uniquey
        push!(unique_indices, findfirst(e -> e == y, envelopey))
    end

    # Filter x and y arrays based on these indices
    envelopex = envelopex[unique_indices]
    envelopey = envelopey[unique_indices]

    extrema_indices = LocalExtremaIndices(envelopey)
    
    println("Local extrema indices: ", extrema_indices)
    println("Local extrema x: ", envelopex[extrema_indices])
    println("Local extrema y: ", envelopey[extrema_indices])

    sum_peaks = Array{Float64}(undef, Int(2 * rabii.j + 1))

    k = 1
    for i in 3:length(extrema_indices)
        if isodd(i)
            sum_peaks[k] = sum(envelopey[findall(e -> envelopex[extrema_indices[i - 2]] < e < envelopex[extrema_indices[i]], envelopex)])
            println("Peak $(k): ", sum_peaks[k])
            k += 1
        end
    end

    println("Peak check: ", sum(sum_peaks))

    if showgraph
        p = scatter(sf[1], log10.(sf[2]),  ylims=(log10(threshold) - 5, 0), xlims=(0,3), title="$(k)", markeralpha=0.5, markerstrokewidth=0, xlabel="\$E\$", ylabel=raw"$\log_{10}S$", label=raw"$|\psi_{\mathrm{GS}}\rangle$")
        p = scatter!(p, envelopex, log10.(envelopey), color=:red, markerstrokewidth=0)
        display(plot(p))
    end
    if savedata
        Export("$(PATH)sf_$(String(rabii))_$(λf)", sf[1], sf[2])
        Export("$(PATH)sf_$(String(rabii))_$(λf)_envelope", envelopex, envelopey)
    end

    return sum_peaks
end
