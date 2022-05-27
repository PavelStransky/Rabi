include("Rabi.jl")

rabi = Rabi(N=500, R=100, λ=0.75)
energies, vectors = eigenstates(rabi)
a, b = ProjectParity(rabi, vectors[1], vectors[2])
pa = projector(a)

rabi = Rabi(rabi; λ=rabi.λ / 2)
operators = AllOperators(rabi)
operators[4] = :P => pa
operators[8] = :r_r => pa
ev = ExpectationValues(rabi, operators; Ψ0=a, saveData = false)
