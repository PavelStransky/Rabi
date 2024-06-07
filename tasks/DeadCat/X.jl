using Distributed

@everywhere const WORKERS = 4

include("../../Calculation.jl")
pyplot(size=(1200, 1200))
# const PATH = "d:/results/Rabi/deadcat/f3/"
const PATH = ""

R = 10000.0
μ = 0.132 / R
# μ = 1.43 / R

rabi = Rabi(R=R, λ=0.75, δ=0.5, μ=μ, ν=μ)
result = ExpectationValues(rabi, [:x => X(rabi)]; mint=0.0, maxt=100.0, numt=1001, showGraph=true, saveGraph=false, saveData=true, asymptotics=false)
