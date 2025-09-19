include("../../Rabi.jl")

const PATH = "d:/results/rabi/Schnellbruder/"

pyplot(size=(1000, 800))

R = 20
λi = 1.5
λf = -sqrt(2.0) / 5
# λf = -sqrt(209) / 12
λf = -1.5

δ = 0.5

n = 300
j = 1 // 2

# Initial and final systems
rabii = Rabi(N=n, R=R, λ=λi, δ=δ, j=j)
rabif = Rabi(N=n, R=R, λ=λf, δ=δ, j=j)

gs = SingleWellState(rabii)	

limits = 1.5
wignerMesh = 151

clim = (-0.08, 0.08)
firstIndex = 1

Wigner(rabif; Ψ0=gs, operators=[:Jx=>Jx(rabif), :q=>X(rabif), :Jy=>Jy(rabif), :p=>P(rabif), :Jz=>Jz(rabif), :P=>projector(gs)], operatorsColor=[:gray, :red, :gray, :red, :gray, :blue], operatorsMarker=[:diamond, :square, :diamond, :square, :diamond, :pentagon],
operatorsLayout=(7, 2), firstIndex=firstIndex, lastIndex=-1,
maxt=6*3.81014, numt=12, xs=LinRange(-limits, limits, wignerMesh), ys=LinRange(-limits, limits, wignerMesh), marginals=true,
clim=clim, saveData=true, saveGraph=true, showGraph=true)