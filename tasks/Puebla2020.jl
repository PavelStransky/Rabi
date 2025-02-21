include("../Rabi.jl")

const PATH = "d:/results/rabi/wf/"

pyplot(size=(600, 500))

R = 20
λi = 1.0
λf = 0.25

showGraph = false
wignerMesh = 401

# Level Dynamics
_, pld = LevelDynamics(Rabi(R=R,N=150), ps=LinRange(0, 1, 301), limit=150, ylims=(-1,2), saveGraph=true)

# Initial and final systems
rabii = Rabi(R=R, λ=λi)
rabif = Copy(rabii; λ=λf)

ei, vi = eigenstates(rabii)
gs = vi[1]

println(ExpectationValue("E", H(rabif), gs, rabif) / R)

# Exact parity states
a, b = ProjectParity(rabii, vi[1], vi[2])
c = sqrt(0.5) * (a + b)
d = sqrt(0.5) * (a - b)

# Asymptotic parity state (Puebla2020)
as = ΨGSP(rabii)

### All states are normalized

# Strength function
sfgs = StrengthFunction(rabif, gs) 
sfa = StrengthFunction(rabif, a)
sfc = StrengthFunction(rabif, c)
sfas = StrengthFunction(rabif, as)

""" 
    Positive parity states
"""
function positiveParity(sf)
    e, s, p = sf
    return e[p .> 0], s[p .> 0]
end

""" 
    Negative parity states
"""
function negativeParity(sf)
    e, s, p = sf
    return e[p .< 0], s[p .< 0]
end


# p = scatter(positiveParity(sfgs), markeralpha=0.5, markerstrokewidth=0, xlabel="\$E\$", ylabel="\$S\$", label=raw"$|\psi_{\mathrm{GS}}\rangle$")
# p = scatter!(p, sfc[1:2], markeralpha=0.5, markerstrokewidth=0, label=raw"$|\psi_{\mathrm{GS}}^{+}\rangle$")
# p = scatter!(p, sfas[1:2], smarkeralpha=0.5, markerstrokewidth=0, label=raw"$|\psi_{\infty}^{+}\rangle$")

p = scatter(positiveParity(sfgs), xlims=(sfgs[1][1], 3), ylims=(-0.001, 0.1), title="Positive parity", markeralpha=0.5, markerstrokewidth=0, xlabel="\$E\$", ylabel="\$S\$", label=raw"$|\psi_{\mathrm{GS}}\rangle$")
p = scatter!(p, positiveParity(sfa), markeralpha=0.5, markerstrokewidth=0, label=raw"$|\psi_{\mathrm{GS}}^{+}\rangle$")
p = scatter!(p, positiveParity(sfc), markeralpha=0.5, markerstrokewidth=0, label=raw"$|\psi_{\mathrm{GS}}^{R}\rangle$")
p = scatter!(p, positiveParity(sfas), smarkeralpha=0.5, markerstrokewidth=0, label=raw"$|\psi_{\infty}^{R}\rangle$")
display(plot(p))

savefig(p, "$(PATH)sf_positive($(rabif.R),$(rabif.λ),$(rabif.δ),$(rabif.μ)).png")

p = scatter(negativeParity(sfgs), xlims=(sfgs[1][1], 3), ylims=(-0.001, 0.1), title="Negative parity", markeralpha=0.5, markerstrokewidth=0, xlabel="\$E\$", ylabel="\$S\$", label=raw"$|\psi_{\mathrm{GS}}\rangle$")
p = scatter!(p, negativeParity(sfa), markeralpha=0.5, markerstrokewidth=0, label=raw"$|\psi_{\mathrm{GS}}^{+}\rangle$")
p = scatter!(p, negativeParity(sfc), markeralpha=0.5, markerstrokewidth=0, label=raw"$|\psi_{\mathrm{GS}}^{R}\rangle$")
p = scatter!(p, negativeParity(sfas), smarkeralpha=0.5, markerstrokewidth=0, label=raw"$|\psi_{\infty}^{R}\rangle$")
display(plot(p))

savefig(p, "$(PATH)sf_negative($(rabif.R),$(rabif.λ),$(rabif.δ),$(rabif.μ)).png")

e, sf = sfas
sf = λf .- 2.5 .* sf
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
#Wigner(rabif; Ψ0=a, operators=[:r_r=>pa], ts=LinRange(0, 20, 1001), index=1:50:1001, xs=LinRange(-1, 1, wignerMesh), ys=LinRange(-1, 1, wignerMesh), clim=(-0.1, 0.1), saveData=false, showGraph=showGraph)
Wigner(rabif; Ψ0=as, operators=[:r_r=>pas], ts=LinRange(0, 20, 1001), index=1:10:1001, xs=LinRange(-1, 1, wignerMesh), ys=LinRange(-1, 1, wignerMesh), clim=(-0.15, 0.15), saveData=false, showGraph=showGraph)

