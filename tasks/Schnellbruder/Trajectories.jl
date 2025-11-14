using DifferentialEquations

include("../../Rabi.jl")

const PATH = "d:/results/rabi/Schnellbruder/"

pyplot(size=(1000, 800))

timespan = pi * 2 / 3
timespan = 3 * pi

R = 20
λi = 1.5
λf = -sqrt(2.0) / 5
# λf = -sqrt(209) / 12
# λf = -1.5

δ = 0.5

n = 300
j = 4 // 2

# Initial and final systems
rabii = Rabi(N=n, R=R, λ=λi, δ=δ, j=j)
rabif = Rabi(N=n, R=R, λ=λf, δ=δ, j=j)

function EquationOfMotion!(dx, x, parameters, t)
    q, p = x
    rabi, m = parameters

    s2 = 2 * rabi.λ^2 * (q^2 + rabi.δ^2 * p^2) + 1    
    s = sqrt(s2)

    dx[1] = p * (1 + rabi.λ^2 * rabi.δ^2 * m / rabi.j / s)
    dx[2] = -q * (1 + rabi.λ^2 * m / rabi.j / s)
end

function Trajectory(rabif, t, λi)
    x0 = [-sqrt(0.5 * (λi^2 - 1 / λi^2)), 0]
    
    timeInterval = (0.0, t)
    solver = TsitPap8()
    tolerance = 1E-10
    fnc = ODEFunction(EquationOfMotion!)

    colours = palette(:auto)
    colourIndex = 3

    p = plot()

    for m = -rabif.j:rabif.j
        problem = ODEProblem(fnc, x0, timeInterval, (rabif, m))

        saveat = collect(range(max(t - timespan, 0), t, step=0.01))
        time = @elapsed solution = solve(problem, solver, reltol=tolerance, abstol=tolerance, verbose=true, saveat=saveat)

        println("m = $m, time = $(time)s, solution = $(length(solution))")

        p = plot!(p, solution, idxs=(1,2), lw=2, alpha=0.4, c=colours[colourIndex], label="m=$m", xlabel=raw"$q$", ylabel=raw"$p$")
        p = scatter!(p, [solution(t)[1]], [solution(t)[2]], markersize=10, c=colours[colourIndex], label=nothing)

        colourIndex += 1

        Export(PATH * "trajectory_rabi_$(rabif.λ)_$(t)_$(Int[2*m])", solution[1,:], solution[2,:])
    end

    display(p)
end

Trajectory(rabif, pi * 8 / 3, λi)