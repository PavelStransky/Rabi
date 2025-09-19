include("../../Rabi.jl")

const PATH = "d:/results/rabi/Schnellbruder/"

pyplot(size=(1000, 800))

R = 20
λi = 1.5
λf = -sqrt(2.0) / 5
δ = 0.5

n = 250
i = 16

# Initial and final systems
rabii = Rabi(N=n, R=R, λ=λi, δ=δ)
rabif = Rabi(N=n, R=R, λ=λf, δ=δ)

λs = LinRange(-0.3, 2.75, 751)
λs = LinRange(-2.5, -0.25, 555)

# Level Dynamics
ld, pld = LevelDynamics(rabii, ps=λs, limit=200, ylims=(-1,2.0))
Export("$(PATH)ld_($(rabii.R),$(rabii.δ))", λs, ld)

energies, vectors = eigenstates(rabii)
Ψ = SingleWellState(rabii)

print(ExpectationValue(H(rabif), Ψ) / rabif.R)

# Strength function
sf = StrengthFunction(rabif, Ψ)
ps = StrengthFunction(rabii; λf=λf, showgraph=true, envelope_window=Int(6 * rabii.j + 2), savedata=true)

p = scatter(sf[1], sf[2], markeralpha=0.7, markerstrokewidth=0, markersize=6, xlims=(-2, 2))
display(p)

Export("$(PATH)sf_($(rabii.R),$(rabii.λ),$(rabif.λ),$(rabii.δ))", sf[1], sf[2])

_, p = Overlap(vectors; limit=length(vectors) ÷ 3)
savefig(p, "$(PATH)psi_($(rabii.R),$(rabii.δ)).png")

p = scatter(pld, λf .- 1.0 .* sf[2], sf[1], markeralpha=0.9, markerstrokewidth=0, markersize=8)
p = plot(p, [λi, λf], [energies[1], ExpectationValue(H(rabif), Ψ) / rabif.R], linewidth=3, arrow=(:closed, 1.0))
display(plot(p))

savefig(p, "$(PATH)sf_($(rabii.R),$(rabii.λ),$(rabif.λ),$(rabii.δ)).png")