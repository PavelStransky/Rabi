include("../../Export.jl")
include("../../Rabi.jl")

pyplot(size=(1000, 1000))

BLAS.set_num_threads(1)     # To prevent the Stack Overflow error

const PATH = "d:/results/Rabi/deadcat/X"

R = 100.0
μ = 0.13 / R
ν = μ

c = 0.5
d = 0.02

mint = 0.0
maxt = 500.0
numt = 1000

rabi = Rabi(R=R, λ=0.75, δ=0.5, μ=μ, ν=ν)
result = ExpectationValues(rabi, [:Jx=>Jx(rabi), :Jy=>Jy(rabi), :Jz=>Jz(rabi), :q=>X(rabi), :p=>P(rabi), :n=>N(rabi)]; mint=mint, maxt=maxt, numt=numt, showGraph=true, saveGraph=true, saveData=true, asymptotics=false)
result = ExpectationValuesLindblad(rabi, [:Jx=>Jx(rabi), :Jy=>Jy(rabi), :Jz=>Jz(rabi), :q=>X(rabi), :p=>P(rabi), :n=>N(rabi)], [c * Jm(rabi)]; mint=mint, maxt=maxt, numt=numt, showGraph=true, saveGraph=true, fname="Jm_", saveData=true, asymptotics=false)
result = ExpectationValuesLindblad(rabi, [:Jx=>Jx(rabi), :Jy=>Jy(rabi), :Jz=>Jz(rabi), :q=>X(rabi), :p=>P(rabi), :n=>N(rabi)], [c * Jp(rabi)]; mint=mint, maxt=maxt, numt=numt, showGraph=true, saveGraph=true, fname="Jp_", saveData=true, asymptotics=false)
result = ExpectationValuesLindblad(rabi, [:Jx=>Jx(rabi), :Jy=>Jy(rabi), :Jz=>Jz(rabi), :q=>X(rabi), :p=>P(rabi), :n=>N(rabi)], [d * A(rabi)]; mint=mint, maxt=maxt, numt=numt, showGraph=true, saveGraph=true, fname="A_", saveData=true, asymptotics=false)
result = ExpectationValuesLindblad(rabi, [:Jx=>Jx(rabi), :Jy=>Jy(rabi), :Jz=>Jz(rabi), :q=>X(rabi), :p=>P(rabi), :n=>N(rabi)], [d * Ad(rabi)]; mint=mint, maxt=maxt, numt=numt, showGraph=true, saveGraph=true, fname="Ad_", saveData=true, asymptotics=false)
