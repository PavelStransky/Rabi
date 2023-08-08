const PATH = "d:\\results\\Lipkin\\"

include("../Lipkin.jl")

pyplot(size=(1200, 1000))

# Habilitace
lipkin = Lipkin(j=10)
x = LinRange(0, 1, 1000)
e, p = LevelDynamics(lipkin; ps=x, type = :ξ, saveGraph=true)
Export("d:\\qpt2", x, e)

lipkin = Lipkin(j=10, ξ=0.6)
x = LinRange(-5, 5, 2000)
e, p = LevelDynamics(lipkin; ps=x, type = :ζ, saveGraph=true, ylims=(-2, 1))
Export("d:\\qpt1", x, e)

lipkin = Lipkin(j=30, ξ=0.6)
x = LinRange(-3, 3, 2000)
e, p = LevelDynamics(lipkin; ps=x, type = :ζ, saveGraph=true, ylims=(-2, 1))
