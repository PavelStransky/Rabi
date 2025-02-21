const PATH = ""

include("../SNAIL.jl")
include("../Export.jl")

BLAS.set_num_threads(1)     # To prevent the Stack Overflow error

pyplot(size=(1200, 1000))

snail = SNAIL(Δ=0, K=4, N=100)
x = LinRange(0, 25, 100)
e, p = LevelDynamics(snail; ps=x, type = :ϵ2, saveGraph=false, ylims=(-600, 600))

snail = SNAIL(Δ=0, K=4, ϵ2=15, N=100)
result = ExpectationValuesLindblad(snail, [:n=>N(snail)], [0.1 * A(snail)]; mint=0, maxt=100, numt=100, showGraph=true, saveGraph=true, fname="Jm_", saveData=true, asymptotics=false, solver=TsitPap8(), maxiters=10000)
