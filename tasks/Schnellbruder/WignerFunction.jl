include("../../Rabi.jl")

const PATH = "d:/results/rabi/Schnellbruder/"

gr()
default(size=(1920,1080))

R = 50
λi = 0
# λf = -sqrt(2.0) / 5
# λf = -sqrt(209) / 12
# λf = -1.5
λf = 1.5
μ = 0.13 / R

δ = 0.5

n = 300
j = 1 // 2

# Initial and final systems
rabii = Rabi(N=n, R=R, λ=λi, δ=δ, j=j, μ=μ, ν=μ)
rabif = Rabi(N=n, R=R, λ=λf, δ=δ, j=j, μ=μ, ν=μ)

gs = ΨGS(rabii)	

limits = 1.5
wignerMesh = 501

clim = (-0.08, 0.08)
firstIndex = 1

# Wigner(rabif; Ψ0=gs, operators=[:Jx=>Jx(rabif), :q=>X(rabif), :Jy=>Jy(rabif), :p=>P(rabif), :Jz=>Jz(rabif), :P=>projector(gs)], operatorsColor=[:gray, :red, :gray, :red, :gray, :blue], operatorsMarker=[:diamond, :square, :diamond, :square, :diamond, :pentagon],
# operatorsLayout=(7, 2), firstIndex=firstIndex, lastIndex=-1,
# maxt=20*pi, numt=60, xs=LinRange(-limits, limits, wignerMesh), ys=LinRange(-limits, limits, wignerMesh), marginals=false,
# clim=clim, saveData=true, saveGraph=true, showGraph=false)

firstIndex = 1
clim = (-0.15, 0.15)
Wigner(rabif; Ψ0=gs, operators=[], firstIndex=firstIndex, lastIndex=-1,
maxt=60, numt=3000, xs=LinRange(-2, 2, wignerMesh), ys=LinRange(-1, 1, wignerMesh), marginals=true,
clim=clim, saveData=false, saveGraph=true, showGraph=false)