using Distributed
@everywhere const WORKERS = 16
include("../../Calculation.jl")
@everywhere const PATH = "d:/results/Rabi/deadcat/minmax/"
# @everywhere const PATH = ""

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

function Calculate(showGraph=false)
    Rs = [50 100 200 300 400 500 800 1000]
    Rs = [300 400 500 800 1000]
    Rs = [100]

    λ = 0.55
    λ = 0.75

    mint = 20.0
    maxt = 30.0
    numt = 4000

    for R in Rs
        minx = 0.0
        maxx = 12 * 0.57 / R
        maxx = 12 * 0.13 / R
        numx = 3000
        μs = LinRange(minx, maxx, numx + 1)
    
        rabi = Rabi(R=R, λ=λ, δ=0.5)

        result = DQPT(rabi, μs=μs, mint=mint, maxt=maxt, numt=numt + 1, showGraph=false, parallel=true, saveData=false)
        mins, maxs = FindExtremes(result, 5, margin=maximum([trunc(Int, 2000 / R), 20]))

        xs = [(maxt - mint) / numt * mins[j][2] + mint for j = 1:length(mins)]
        ys = [(maxx - minx) / numx * mins[j][1] + minx for j = 1:length(mins)]
        zs = [mins[j][3] for j = 1:length(mins)]

        if showGraph
            p = scatter(xs, ys, title="R=$(R)", legend=false)
            display(p)
        end

        Export("$(PATH)lambda=$(λ)_mins_R=$(R).txt", xs, ys, zs)

        xs = [(maxt - mint) / numt * maxs[j][2] + mint for j = 1:length(maxs)]
        ys = [(maxx - minx) / numx * maxs[j][1] + minx for j = 1:length(maxs)]
        zs = [maxs[j][3] for j = 1:length(maxs)]

        if showGraph
            p = scatter(xs, ys, title="R=$(R)", legend=false)
            display(p)
        end

        Export("$(PATH)lambda=$(λ)_maxs_R=$(R).txt", xs, ys, zs)
    end
end

Calculate()