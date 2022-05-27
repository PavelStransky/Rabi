using Distributed

@everywhere const WORKERS = 8

include("../Calculation.jl")

@everywhere const PATH = ""

pyplot(size=(1000, 800))

rabi = Rabi(R=100, λ=0.75, δ=0.5, μ=0.001, ν=0.001)
ExpectationValues(rabi, AllOperators(rabi); saveData=false)
c = Compute(rabi; λs=LinRange(0, 1.5, 601), parallel=true)

rabi = Rabi(rabi; ν=0.0)
ExpectationValues(rabi, AllOperators(rabi); saveData=false)
c = Compute(rabi; λs=LinRange(0, 1.5, 601), parallel=true)

rabi = Rabi(R=100, λ=0.75, δ=0.5, μ=0.4, ν=0.4)
ExpectationValues(rabi, AllOperators(rabi); saveData=false)
c = Compute(rabi; λs=LinRange(0, 1.5, 601), parallel=true)

rabi = Rabi(rabi; ν=0.0)
ExpectationValues(rabi, AllOperators(rabi); saveData=false, maxt=200)
c = Compute(rabi; λs=LinRange(0, 1.5, 601), parallel=true)
