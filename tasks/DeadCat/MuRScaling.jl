using Distributed
@everywhere const WORKERS = 16
include("../../Calculation.jl")
@everywhere const PATH = "d:/results/Rabi/deadcat/"

minx = 0.0
maxx = 0.05
numx = 1000
μs = LinRange(minx, maxx, numx + 1)

mint = 15.0
maxt = 25.0
numt = 1000

function FindExtremes(result, k; margin=1)
    mins = []

    for i in (margin + 1):(size(result, 2) - margin), j in (margin + 1):(size(result, 3) - margin)
        minimum = true
        for m = -margin:margin, n = -margin:margin
            if result[k, i, j] > result[k, i + m, j + n]
                minimum = false
                break
            end
        end
        if minimum
            push!(mins, [i, j, result[k, i, j]])
        end
    end
    return mins
end

function Export(fname, xs, ys, zs)
    open(fname, "w") do io
        for (x, y, z) in zip(xs, ys, zs)
            println(io, "$(x)\t$(y)\t$(z)")
        end
    end
end

function Calculate()
    R = 25

    while R <= 500
        rabi = Rabi(R=R, λ=0.75, δ=0.5)

        result = DQPT(rabi, μs=μs, mint=mint, maxt=maxt, numt=numt + 1, showGraph=false, parallel=true, saveData=false)
        mins = FindExtremes(result, 5, margin=trunc(Int, 1000 / R))

        xs = [(maxt - mint) / numt * mins[j][2] + mint for j = 1:length(mins)]
        ys = [(maxx - minx) / numx * mins[j][1] + minx for j = 1:length(mins)]
        zs = [mins[j][3] for j = 1:length(mins)]

        p = scatter(xs, ys, title="R=$(R)", legend=false)
        display(p)

        Export("$(PATH)extremes_$R.txt", xs, ys, zs)

        R *= 2
    end
end

Calculate()