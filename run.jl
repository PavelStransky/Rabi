include("Export.jl")
include("Rabi.jl")

pyplot(size=(1000, 1000))

BLAS.set_num_threads(1)     # To prevent the Stack Overflow error

const PATH = "d:\\results\\Rabi\\"

R = 3200.0
μ = 0.13 / R
ν = μ

rabi = Rabi(R=R, λ=0.75, δ=0.5, μ=μ, ν=ν)
result = ExpectationValues(rabi, AllOperators(rabi); mint=0.0, maxt=50.0, numt=1000, showGraph=true, saveGraph=true, saveData=true, asymptotics=false)

# xs = log10.(μs)
# ts = LinRange(0, 50, 1001)

# p1 = heatmap(ts, xs, result[1,:,:], color=:coolwarm, clims=(-0.1, 0.1), title="Jx")
# p5 = heatmap(ts, xs, result[5,:,:], color=:coolwarm, clims=(-0.5, 0.5), title="q")
# p = plot(p1, p5, layout=(2, 1))
# savefig(p, PATH * "jx q $R.png")

# p2 = heatmap(ts, xs, result[2,:,:], color=:coolwarm, clims=(-0.1, 0.1), title="Jy")
# p6 = heatmap(ts, xs, result[6,:,:], color=:coolwarm, clims=(-0.5, 0.5), title="p")
# p = plot(p2, p6, layout=(2, 1))
# savefig(p, PATH * "jy p $R.png")

# p3 = heatmap(ts, xs, result[3,:,:], clims=(-0.5, 0), title="Jz")
# p4 = heatmap(ts, xs, result[4,:,:], clims=(0, 1), title="Survival probability")
# p = plot(p3, p4, layout=(2, 1))
# savefig(p, PATH * "jz P $R.png")
