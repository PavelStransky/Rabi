include("../Lipkin.jl")

print(dense(Jz(Lipkin(j=10))))
print(eigenstates(Lipkin(j=10, λ=0))[1])
print(Size(Lipkin(j=10)))

lipkin = Lipkin(j=50, λ=1)
LevelDynamics(lipkin; ps=LinRange(0, 3, 100), ylims=(-1.7, -0.6))
LevelDynamics(lipkin; ps=LinRange(0, 1, 100), type = :ω)