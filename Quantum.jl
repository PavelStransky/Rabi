using Plots
using QuantumOptics
using Statistics
using LinearAlgebra
using Formatting

abstract type QuantumSystem end

PGS(system::QuantumSystem) = projector(ΨGS(system))

""" Expectation value of a given operator """
function ExpectationValue(name, operator, Ψ, system::QuantumSystem=nothing)
    name = String(name)
    result = real(expect(operator, Ψ))
    if occursin("_log", name) result = log10.(result)
    elseif occursin("_r", name) result = -log.(abs.(result)) / Size(system) end 
    return result
end

ExpectationValue(operator, Ψ) = real(expect(operator, Ψ))

function Label(name) 
    name = String(name)

    if occursin("_log", name)
        result = "\$\\log" * replace(name, "_log"=>"") * "\$"
    elseif occursin("_r", name)
        result = "\$" * replace(name, "_r"=>"") * "\$"
    else
        result = "\$$name\$"
    end
end

function StrengthFunction(rabif, keti)
    ef, vf = eigenstates(rabif)
    sf = [abs(dagger(vf[j]) * keti)  for (j, _) in enumerate(vf)]
    parity = [ExpectationValue(Parity(rabif), vf[j]) for (j, _) in enumerate(vf)]

    return ef, sf, parity
end

function Overlap(vectors; limit=nothing)
    size = limit === nothing ? length(vectors) : limit

    result = Matrix{Float64}(undef, size, size)
    va = [abs.(vectors[j].data) for j in 1:size]

    for i in 1:size
        for j in i:size
            result[i, j] = sum(va[i] .* va[j])
            result[j, i] = result[i, j]
        end
    end

    p = heatmap(result, color=:nipy_spectral, clim=(0, 1))
    display(p)

    return result, p
end

" Wigner function "
function Wigner(system::QuantumSystem; husimi=false, Ψ0=nothing, ts=[0,1,2], index=nothing, operators=[], xs=LinRange(-1, 1, 101), ys=nothing, showGraph=true, saveData=true, saveGraph=true, log=false, clim=(-0.2, 0.2), kwargs...)
    pyplot(size = (1000, 1000))

    # Range in y direction
    if ys === nothing ys = xs end
    
    # Initial state
    if Ψ0 === nothing Ψ0 = ΨGS(system) end

    # Schroedinger time evolution
    print("Schrödinger ", system, "...")
    time = @elapsed tout, Ψt = timeevolution.schroedinger(ts, Ψ0, H(system))
    println(time, "s")

    opvalues = [name=>[ExpectationValue(name, operator, Ψ, system) for Ψ in Ψt] for (name, operator) in operators]

    xss = xs * sqrt(2.0 * Size(system))
    yss = ys * sqrt(2.0 * Size(system))

    if index === nothing index = 1:length(tout) end
    title = husimi ? "Husimi" : "Wigner"

    for i = index
        tstr = format(tout[i], precision=2)

        print("T[$i]=$tstr trace...")
        time = @elapsed ρ = DensityMatrix(system, Ψt[i])
        print("$time $title...")

        time = @elapsed w = husimi ? real.(qfunc(ρ, xss, yss)) : wigner(ρ, xss, yss)
        print("$time display...")

        logstr = ""
        if log
            w = -sign.(w) .* log10.(abs.(w))
            replace!(x -> (x > clim[1] && x < clim[2]) ? x : 0.0, w)
            replace!(x -> x < 0.0 ? clim[1] - x : x, w)
            replace!(x -> x > 0.0 ? clim[2] - x : x, w)
            logstr = "_log"
        end            

        time = @elapsed begin
            p = heatmap(xs, ys, transpose(w), color=:bwr, grid=false, title="$title $system, \$t=$tstr\$", xlabel=raw"$q$", ylabel=raw"$p$", clim=clim, kwargs...)

            len = length(opvalues)
            psp = Array{Any}(undef, len)
            pspi = 1
            for (name, opvalue) in opvalues
                psp[pspi] = plot(tout, opvalue)
                psp[pspi] = scatter!(psp[pspi], [tout[i]], [opvalue[i]], markersize=10, markeralpha=0.7, xlabel=raw"$t$", ylabel=Label(name), legend=false)
                pspi += 1
            end

            if len > 0
                ps = len == 1 ? plot(psp[1]) : plot(psp..., layout=grid(1, len))
                p = plot(ps, p, layout=grid(2, 1, heights=[0.2 ,0.8]))
            end

            showGraph && display(plot(p))
            saveGraph && savefig(p, "$(PATH)$(title)_$(system)_$i$logstr.png")
        end

        println(time)

        saveData && Export("$(PATH)$(title)_$(system)_$tstr$logstr.txt", xs, ys, w)
    end

    return true
end

""" Expectation values of specific operator """
function ExpectationValues(system::QuantumSystem, operators; Ψ0=nothing, mint=0.0, maxt=100.0, numt=1001, showGraph=true, saveGraph=true, saveData=true, asymptotics=true, kwargs...)    
    # Initial state
    if Ψ0 === nothing Ψ0 = ΨGS(system) end

    print("Expectation values ", system, " E=", real(expect(H(system), Ψ0)) / Size(system), "...")

    numop = length(operators)
    result = Array{Float64}(undef, numop, numt)

    time = @elapsed begin
        u = U(system, mint)
        du = U(system, (maxt - mint) / numt)
        Ψ = u * Ψ0      # Initial state

        for i in 1:numt
            for (j, (name, operator)) in enumerate(operators)
                result[j, i] = ExpectationValue(name, operator, Ψ, system)
            end
            Ψ = du * Ψ
        end
    end
    println(time)

    if asymptotics
        asymptoticValues = AsymptoticValues(system, operators; Ψ0=Ψ0)
    end

    tout = LinRange(mint, maxt, numt)

    if showGraph || saveGraph
        ps = Array{Any}(undef, numop)
        for (j, (name, _)) in enumerate(operators)
            ps[j] = plot(tout, result[j,:], xlabel="\$t\$", title=Label(name), legend=false)
            if asymptotics
                ps[j] = plot!(ps[j], [mint, maxt], [asymptoticValues[j], asymptoticValues[j]], color = :red, width=2)
            end
        end

        p = plot(ps..., layout=(2, numop ÷ 2); kwargs...)
        showGraph && display(p)
        saveGraph && savefig(p, "$(PATH)expectation_$system.png")
    end

    if saveData
        for (j, (name, _)) in enumerate(operators)
            Export("$(PATH)$(name)_$system.txt", tout, result[j,:], asymptotics ? asymptoticValues[j] : nothing)
        end
    end

    return result
end

""" Asymptotic values of the operators """
function AsymptoticValues(system::QuantumSystem, operators; Ψ0=nothing, mint=200.0, maxt=1000.0, numt=2000)
    print("Asymptotic ")

    meant = ExpectationValues(system, operators; Ψ0=Ψ0, mint=mint, maxt=maxt, numt=numt, showGraph=false, saveGraph=false, saveData=false, asymptotics=false)

    asymptotic = [mean(meant[j,:]) for (j, _) in enumerate(operators)]

    #= Coherence and Purity
    if length(operators) >= 3
        print("CP...")
        append!(asymptotic, mean(meant[1,:].^2 .+ meant[2,:].^2))
        append!(asymptotic, 0.5 + 2 * mean(meant[1,:].^2 .+ meant[2,:].^2 .+ meant[3,:].^2))
    end
    =#
    println(asymptotic)

    return asymptotic
end

""" Asymptotic values of the operators """
function AsymptoticValuesMatrix(system::QuantumSystem, operators; Ψ0=nothing)
    # Initial state
    if Ψ0 === nothing Ψ0 = ΨGS(system) end

    print("$system...")

    result = Array{Float64}(undef, 0)

    time = @elapsed begin
        energies, vectors = eigenstates(system)

        p = [abs(dagger(vector) * Ψ0)^2 for vector in vectors]

        for (name, operator) in operators
            a = [real(dagger(vector) * operator * vector) for vector in vectors]
            append!(result, sum(a .* p))
        end
    end
    println(time)

    return result
end

""" Expectation values of various operators """
function SurvivalLog(system::QuantumSystem; mint=1e-1, maxt=100, numt=10000, Ψ0=nothing)    
    # Initial state
    if Ψ0 === nothing Ψ0 = ΨGS(system) end

    PΨ = projector(Ψ0)
    St = Array{Float64}(undef, numt)

    # Schroedinger time evolution
    print("Evolution...")
    time = @elapsed begin
        u = exp(dense(-im * H(system) * mint))
        du = exp(dense(-im * H(system) * (maxt - mint) / numt))
        Ψ = u * Ψ0      # Initial state

        for i in 1:numt
            St[i] = real(expect(PΨ, Ψ))
            Ψ = du * Ψ
        end
    end
    print(time)
    
    tout = LinRange(mint, maxt, numt)
    p = plot(tout, St, title="S log", xaxis=:log, yaxis=:log)
    display(p)    

    Export("d:\\Sl.txt", log10.(tout), log10.(St), 0)
end

""" Level dynamics """
function LevelDynamics(system::QuantumSystem; ps=LinRange(0, 2, 401), type=:λ, limit=200, saveGraph=false, kwargs...)
    print("Level dynamics in $type: ", system, "...")

    limit = min(Dimension(system), limit)

    time = @elapsed begin
        result = Array{Float64}(undef, 0)    
        for p in ps
            system = Copy(system; N=-limit, type=>p)
            energies, _ = eigenstates(system, limit)
            append!(result, energies[1:limit])
        end
    end
    println(time, "s")
    
    result = transpose(reshape(result, limit, :))
    p = plot(ps, result, legend=false, color=RGBA(0, 0, 0, 0.5), xlabel=Label(type), ylabel=raw"$E$"; kwargs...)
    display(p)

    saveGraph && savefig(p, "$(PATH)ld_$system.png")
    
    return result, p
end