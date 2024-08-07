include("../Rabi.jl")

const PATH = ""

pyplot(size=(1000, 800))

R = 30
λ = 0.75
μ = 0.0
ν = μ

n = 250
i = 16

# Initial and final systems
rabi = Rabi(N=n, R=R, λ=λ, μ=μ, ν=ν)

# Level Dynamics
_, pld = LevelDynamics(rabi, ps=LinRange(0, 1.2, 1201), limit=200, ylims=(-1,1.2))

energies, vectors = eigenstates(rabi)
Ψ = ΨGS(rabi)

# Ground state
rabi0 = Rabi(N=n, R=R)
e, Ψ = eigenstates(rabi0)

Ψ = Ψ[i]

# Strength function
sf = StrengthFunction(rabi, Ψ)

p = scatter(sf[1], sf[2], markeralpha=0.7, markerstrokewidth=0, markersize=6, xlims=(-2, 2))
display(p)

_, p = Overlap(vectors; limit=length(vectors) ÷ 3)
savefig(p, "$(PATH)psi_($(rabi.R),$(rabi.μ)).png")

p = scatter(pld, λ .- 1.0 .* sf[2], sf[1], markeralpha=0.9, markerstrokewidth=0, markersize=8)
p = plot(p, [0, λ], [e[i], ExpectationValue(H(rabi), Ψ) / rabi.R], linewidth=3, arrow=(:closed, 1.0))
display(plot(p))

savefig(p, "$(PATH)sf_($(rabi.R),$(rabi.λ),$(rabi.δ),$(rabi.μ)).png")

