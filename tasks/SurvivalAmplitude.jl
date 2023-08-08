using Distributed

@everywhere const WORKERS = 8

include("../Calculation.jl")

@everywhere const PATH = "d:\\"

pyplot(size=(1000, 800))

rabi = Rabi(N=100, R=8, λ=0.75, δ=0.5)
SurvivalAmplitude(rabi, maxt=20, numt=10001)