include("../../Rabi.jl")

const PATH = "d:/results/Rabi/deadcat/"

function Compute()
    result = []

    for d = ds
        # ev = ExpectationValues(rabi, [:q=>X(rabi)]; mint=mint, maxt=maxt, numt=numt, showGraph=false, saveGraph=false, fname="A_", saveData=false, asymptotics=false)
        ev = ExpectationValuesLindblad(rabi, [:q=>X(rabi)], [d * A(rabi)]; mint=mint, maxt=maxt, numt=numt, showGraph=true, saveGraph=false, fname="A_", saveData=true, asymptotics=false)
        ev = [ev[j] for j = 1:length(ev)]
        (value, index) = findmin(ev)

        display(plot(ev))

        println("Minimum d = $d, value = $value")

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

R = 100.0
μ = 0.135341379 / R
ν = μ

ds = LinRange(0, 0.1, 21)

mint = 0.0
maxt = 25.0
numt = 1000

rabi = Rabi(R=R, λ=0.75, δ=0.5, μ=μ, ν=ν)

result = Compute()
Save("$(PATH)lindblad_$(R).txt", result)