using Distributed
@everywhere const WORKERS = 16
include("../../Calculation.jl")
@everywhere const PATH = "d:/results/Rabi/deadcat/"

function FindExtremes(result, k; margin=1)
    mins = []
    maxs = []

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

        maximum = true
        for m = -margin:margin, n = -margin:margin
            if result[k, i, j] < result[k, i + m, j + n]
                maximum = false
                break
            end
        end
        if maximum
            push!(maxs, [i, j, result[k, i, j]])
        end
    end

    return mins, maxs
end

function Export(fname, xs, ys, zs)
    open(fname, "w") do io
        for (x, y, z) in zip(xs, ys, zs)
            println(io, "$(x)\t$(y)\t$(z)")
        end
    end
end

function Calculate()
    Rs = [50 100 200 500 1000 2000 5000]

    for R in Rs
        minx = 0.0
        maxx = 10 * 0.135 / R
        numx = 1000
        μs = LinRange(minx, maxx, numx + 1)
    
        mint = 16.0
        maxt = 26.0
        numt = 1000

        rabi = Rabi(R=R, λ=0.75, δ=0.5)

        result = DQPT(rabi, μs=μs, mint=mint, maxt=maxt, numt=numt + 1, showGraph=false, parallel=true, saveData=true)
        mins, maxs = FindExtremes(result, 5, margin=maximum([trunc(Int, 1000 / R), 10]))

        xs = [(maxt - mint) / numt * mins[j][2] + mint for j = 1:length(mins)]
        ys = [(maxx - minx) / numx * mins[j][1] + minx for j = 1:length(mins)]
        zs = [mins[j][3] for j = 1:length(mins)]

        p = scatter(xs, ys, title="R=$(R)", legend=false)
        display(p)

        Export("$(PATH)mins_$R.txt", xs, ys, zs)

        xs = [(maxt - mint) / numt * maxs[j][2] + mint for j = 1:length(maxs)]
        ys = [(maxx - minx) / numx * maxs[j][1] + minx for j = 1:length(maxs)]
        zs = [maxs[j][3] for j = 1:length(maxs)]

        p = scatter(xs, ys, title="R=$(R)", legend=false)
        display(p)

        Export("$(PATH)maxs_$R.txt", xs, ys, zs)
    end
end

Calculate()