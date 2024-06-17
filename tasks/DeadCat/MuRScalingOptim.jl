using Optim
const WORKERS = 1
include("../../Calculation.jl")
const PATH = "d:/results/Rabi/deadcat/minmax/"

function Export(fname, xs, ys, zs)
    open(fname, "w") do io
        for (x, y, z) in zip(xs, ys, zs)
            println(io, "$(x)\t$(y)\t$(z)")
        end
    end
end

function Calculate()
    R = 500
    λ = 0.75

    rabi = Rabi(R=R, λ=λ, δ=0.5)
    Ψ0 = ΨGS(rabi)
    operator = X(rabi)

    maxμ = 0.26 / R
    i = 0

    function f(x)
        μ = x[1]
        t = x[2]

        i += 1

        if μ > maxμ || μ < 0 || t < 17 || t > 27
            return 0
        end

        system = Copy(rabi; μ=μ, ν=μ)

        u = U(system, t)
        Ψ = u * Ψ0      # Initial state

        result  = real(expect(operator, Ψ))
        println("$i $(system) t = $t, min = $result")
        return result
    end

    time = @elapsed result = optimize(f, [0.13 / R, 20])
    # result = optimize(f, [0, maxμ], [15, 25], [0.0013, 20.3])

    println("Time: $time")
    println("Minimum: $(result.minimum) at $(result.minimizer)")

    return result
end

result = Calculate()