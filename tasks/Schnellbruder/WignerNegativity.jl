include("../../Rabi.jl")

const PATH = "d:/results/rabi/Schnellbruder/"

R = 20
λi = 1.5
λf = -sqrt(2.0) / 5
# λf = -sqrt(209) / 12
# λf = -1.5

δ = 0.5

n = 300
j = 1 // 2

# Initial and final systems
rabii = Rabi(N=n, R=R, λ=λi, δ=δ, j=j)
rabif = Rabi(N=n, R=R, λ=λf, δ=δ, j=j)

gs = SingleWellState(rabii)	

limits = 2.0
wignerMesh = 400

firstIndex = 1

result = Wigner(rabif; Ψ0=gs, firstIndex=firstIndex, lastIndex=-1,
maxt=100, numt=1000, xs=LinRange(-limits, limits, wignerMesh), ys=LinRange(-limits, limits, wignerMesh), marginals=false,
saveData=true, saveGraph=false, showGraph=false)

ev = ExpectationValues(rabif, [:Jx=>Jx(rabif), :Jy=>Jy(rabif), :Jz=>Jz(rabif)]; Ψ0=gs, numt=2000, asymptotics=false)