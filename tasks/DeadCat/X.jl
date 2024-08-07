include("../../Rabi.jl")
pyplot(size=(1200, 1200))
const PATH = "d:/results/Rabi/deadcat/f3/"
# const PATH = ""

R = 200.0
μ = 0.129888 / R
# μ = 1.43 / R

rabi = Rabi(R=R, λ=0.75, δ=0.5, μ=μ, ν=μ)
result = ExpectationValues(rabi, [:x => X(rabi)]; mint=0.0, maxt=100.0, numt=2001, showGraph=true, saveGraph=false, saveData=true, asymptotics=false)
