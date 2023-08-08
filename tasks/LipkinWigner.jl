const PATH = "d:\\results\\Lipkin\\"

include("../Lipkin.jl")

pyplot(size=(1200, 1000))


lipkin = Lipkin(j=50, λ=1)
LevelDynamics(lipkin; ps=LinRange(0, 3, 100), ylims=(-1.7, -0.6))
LevelDynamics(lipkin; ps=LinRange(0, 1, 100), type = :ω, ylims=(-1, 0.5), saveGraph=true)