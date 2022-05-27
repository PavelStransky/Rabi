include("../Rabi.jl")

const PATH = ""

pyplot(size=(1000, 800))

R = 10
λ = 0.75
μ = 0.4
ν = μ

# Initial and final systems
rabi = Rabi(R=R, λ=λ, μ=μ, ν=ν)

# Level Dynamics
#_, pld = LevelDynamics(rabi, ps=LinRange(0, 1, 1001), limit=300, ylims=(-1,1))

energies, vectors = eigenstates(rabi)
Ψ = ΨGS(rabi)

# Strength function
sf = StrengthFunction(rabi, Ψ)

# "Odd" and "even" parity states
parity = sf[3]
sfp = Array{Float64}(undef, 0)
ep = Array{Float64}(undef, 0)
sfn = Array{Float64}(undef, 0)
en = Array{Float64}(undef, 0)
for (j, p) in enumerate(parity)
    if p > 0
        append!(ep, sf[1][j])
        append!(sfp, sf[2][j])
    else
        append!(en, sf[1][j])
        append!(sfn, sf[2][j])
    end
end

#_, p = Overlap(vectors; limit=length(vectors) ÷ 3)
#savefig(p, "$(PATH)psi_($(rabi.R),$(rabi.μ)).png")

p = scatter(pld, λ .- 1.0 .* sfp, ep, markeralpha=0.9, markerstrokewidth=0, markersize=8)
p = scatter(p, λ .- 1.0 .* sfn, en, markeralpha=0.9, markerstrokewidth=0, markersize=8)
p = plot(p, [0, λ], [-0.5, ExpectationValue(H(rabi), Ψ) / rabi.R], linewidth=3, arrow=(:closed, 1.0))
display(plot(p))

savefig(p, "$(PATH)sf_($(rabi.R),$(rabi.λ),$(rabi.δ),$(rabi.μ)).png")

pyplot(size=(800, 150))

for i in 1:20
    _, pl = PlotΨ(rabi, vectors[i])
    savefig(pl, "$(PATH)psi$i.png")
end

# Expectation values -  What Puebla2020 shows in Figure 3
pΨ = projector(Ψ)
operators = AllOperators(rabi)
operators[4] = :P => pΨ
operators[8] = :r_r => pΨ

pyplot(size=(1000, 800))
ExpectationValues(rabi, operators; Ψ0=Ψ, maxt=20, saveData=false)

# Wigner functions
#Wigner(rabif; Ψ0=a, operators=[:r_r=>pa], ts=LinRange(0, 20, 1001), index=1:50:1001, xs=LinRange(-1, 1, 201), ys=LinRange(-1, 1, 201), clim=(-0.1, 0.1), saveData=false, showGraph=true)
#Wigner(rabif; Ψ0=as, operators=[:r_r=>pas], ts=LinRange(0, 20, 1001), index=1:50:1001, xs=LinRange(-1, 1, 201), ys=LinRange(-1, 1, 201), clim=(-0.1, 0.1), saveData=false, showGraph=true)
