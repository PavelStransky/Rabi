include("../Rabi.jl")

const PATH = ""

pyplot(size=(1000, 800))

R = 30
λi = 0.75
λf = λi / 2

showGraph = false
wignerMesh = 401

# Level Dynamics
_, pld = LevelDynamics(Rabi(R=R), ps=LinRange(0, 1, 1001), limit=300, ylims=(-1,2), saveGraph=true)

# Initial and final systems
rabii = Rabi(R=R, λ=λi)
rabif = Rabi(rabii; λ=λf)

ei, vi = eigenstates(rabii)

# Exact parity states
a, b = ProjectParity(rabii, vi[1], vi[2])
c = sqrt(0.5) * (a + b)
d = sqrt(0.5) * (a - b)

# Asymptotic parity state
as = ΨGSP(rabii)

# Strength function
sfa = StrengthFunction(rabif, a)
sfc = StrengthFunction(rabif, c)
sfas = StrengthFunction(rabif, as)

p = scatter(sfa, markeralpha=0.5, markerstrokewidth=0, xlabel="\$E\$", ylabel="\$S\$", label=raw"$|\psi_{\mathrm{GS}}\rangle$")
p = scatter!(p, sfc, markeralpha=0.5, markerstrokewidth=0, label=raw"$|\psi_{\mathrm{GS}}^{+}\rangle$")
p = scatter!(p, sfas, smarkeralpha=0.5, markerstrokewidth=0, label=raw"$|\psi_{\infty}^{+}\rangle$")
display(plot(p))

savefig(p, "$(PATH)sf_($(rabif.R),$(rabif.λ),$(rabif.δ),$(rabif.μ)).png")

e, sf = sfas
sf = λf .- 1.0 .* sf
p = scatter(pld, sf, e, markeralpha=0.9, markerstrokewidth=0)
p = plot(p, [λi, λf], [ei[1], ExpectationValue(H(rabif), as) / rabif.R], linewidth=3, arrow=(:closed, 1.0))

display(plot(p))

savefig(p, "$(PATH)ld_($(rabif.R),$(rabif.δ),$(rabif.μ)).png")

# Expectation values -  What Puebla2020 shows in Figure 3
pa = projector(a)
operators = AllOperators(rabif)
operators[4] = :P => pa
operators[8] = :r_r => pa

ExpectationValues(rabif, operators; Ψ0=a, maxt=20, saveData=false)

# Expectation values -  What Puebla2020 should show in Figure 3
pas = projector(as)
operators = AllOperators(rabif)
operators[4] = :P => pas
operators[8] = :r_r => pas

ExpectationValues(rabif, operators; Ψ0=as, maxt=20, saveData=false)

# Wigner functions
Wigner(rabif; Ψ0=a, operators=[:r_r=>pa], ts=LinRange(0, 20, 1001), index=1:50:1001, xs=LinRange(-1, 1, wignerMesh), ys=LinRange(-1, 1, wignerMesh), clim=(-0.1, 0.1), saveData=false, showGraph=showGraph)
Wigner(rabif; Ψ0=as, operators=[:r_r=>pas], ts=LinRange(0, 20, 1001), index=1:50:1001, xs=LinRange(-1, 1, wignerMesh), ys=LinRange(-1, 1, wignerMesh), clim=(-0.1, 0.1), saveData=false, showGraph=showGraph)

