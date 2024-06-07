using Distributed

@everywhere  const WORKERS = 16

if nprocs() <= WORKERS
    addprocs(WORKERS + 1 - nprocs())
end

@everywhere include("../../Rabi.jl")

# @everywhere const PATH = "d:/results/Rabi/deadcat/f5/"
@everywhere const PATH = ""

function Compute()
    R = 100.0
    μ = 0.1325 / R
    ν = μ
    
    ds = LinRange(0, 0.15, 151)
    
    mint = 0.0
    maxt = 30.0
    numt = 1001
    
    rabi = Rabi(R=R, λ=0.75, δ=0.5, μ=μ, ν=ν)
    
    time = @elapsed evs = pmap((d)->ExpectationValuesLindblad(rabi, [:q=>X(rabi)], [d * A(rabi)]; mint=mint, maxt=maxt, numt=numt, showGraph=false, saveGraph=false, fname="A_", saveData=true, asymptotics=false), ds)
    println("Elapsed time: $time")

    result = []
    for (index, d) = enumerate(ds)
        ev = evs[index]	
        ev = [ev[j] for j = 1:length(ev)]
        (value, index) = findmin(ev)

        display(plot(ev))
        push!(result, [d, value])
    end

    return result
end

function Save(fname, data)
    open(fname, "w") do io
        for d in data
            println(io, "$(d[1])\t$(d[2])")
        end
    end
end

result = Compute()
Save("$(PATH)lindblad_100.txt", result)