using Plots
using QuantumOptics
using Statistics
using LinearAlgebra
using Formatting
using Measures

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

function StrengthFunction(system::QuantumSystem, keti)
    ef, vf = eigenstates(system)
    sf = [abs(dagger(vf[j]) * keti)^2  for (j, _) in enumerate(vf)]
    parity = [ExpectationValue(Parity(system), vf[j]) for (j, _) in enumerate(vf)]

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

function Wigner(system::QuantumSystem; husimi=false, Ψ0=nothing, maxt=30.0, numt=31, firstIndex=1, lastIndex=-1, operators=[], operatorsColor=nothing, operatorsMarker=nothing, operatorsLayout=nothing, marginals=false, xs=LinRange(-1, 1, 101), ys=nothing, showGraph=true, saveData=true, saveGraph=true, log=false, clim=(-0.2, 0.2), postProcess=nothing, postProcessParams=nothing, kwargs...)
    """ Wigner function """

    # Range in y direction
    if ys === nothing ys = xs end
    
    # Initial state
    if Ψ0 === nothing Ψ0 = ΨGS(system) end

    if operatorsColor === nothing
        operatorsColor = [:blue for _ in 1:length(operators)]
    end

    if operatorsMarker === nothing
        operatorsMarker = [:circle for _ in 1:length(operators)]
    end

    factor = round(Int, 400 / numt)
    if factor < 1 factor = 1 end
    ts = LinRange(0, maxt, numt * factor + 1)

    firstIndex = (firstIndex - 1) * factor + 1

    println("Wigner $(system), maxt=$maxt, numt=$numt, factor=$factor, firstIndex=$firstIndex, lastIndex=$lastIndex")

    # Schroedinger time evolution
    print("Schrödinger ", system, "...")
    time = @elapsed tout, Ψt = timeevolution.schroedinger(ts, Ψ0, H(system))
    println(time, "s")

    opvalues = [name=>[ExpectationValue(name, operator, Ψ, system) for Ψ in Ψt] for (name, operator) in operators]

    xss = xs * sqrt(Size(system))
    yss = ys * sqrt(Size(system))

    title = husimi ? "Husimi" : "Wigner"

    if lastIndex < 0 || lastIndex > numt * factor
        lastIndex = numt * factor
    end

    for i = firstIndex:factor:lastIndex
        tstr = format(tout[i], precision=2)

        print("T[$i]=$tstr")
        time = @elapsed ρ = DensityMatrix(system, Ψt[i])
        print(" trace...$(time)s")

        time = @elapsed w = husimi ? real.(qfunc(ρ, xss, yss)) : wigner(ρ, xss, yss)
        println(" $title...$(time)s")

        logstr = ""
        if log
            w = -sign.(w) .* log10.(abs.(w))
            replace!(x -> (x > clim[1] && x < clim[2]) ? x : 0.0, w)
            replace!(x -> x < 0.0 ? clim[1] - x : x, w)
            replace!(x -> x > 0.0 ? clim[2] - x : x, w)
            logstr = "_log"
        end            

        time = @elapsed begin
            p = heatmap(xs, ys, transpose(w), c=:bwr, grid=false, title="$title $system, t=$tstr", xlabel="q", ylabel="p", clim=clim, 
            left_margin=10mm, bottom_margin=5mm, 
            kwargs...)

            len = length(opvalues)
            psp = Array{Any}(undef, len)
            pspi = 1
            for (name, opvalue) in opvalues
                psp[pspi] = plot(tout, opvalue, color=operatorsColor[pspi])
                psp[pspi] = scatter!(psp[pspi], [tout[i]], [opvalue[i]], color=operatorsColor[pspi], m=operatorsMarker[pspi], markersize=7,  markeralpha=0.7, 
                title=Label(name),
                legend=false)
                # psp[pspi] = annotate!(psp[pspi], 100, -1, text("t", :black, 10))
                pspi += 1
            end

            # psp[len] = plot!(psp[len], bottom_margin=5mm, left_margin=10mm)

            if marginals
                marginal_x = vec(sum(w, dims=2))
                marginal_x[isnan.(marginal_x)] .= 0.0
                marginal_x = marginal_x / (sum(marginal_x) * (xs[2] - xs[1]))

                marginal_y = vec(sum(w, dims=1))
                marginal_y[isnan.(marginal_y)] .= 0.0
                marginal_y = marginal_y / (sum(marginal_y) * (ys[2] - ys[1]))

                p = plot!(p, xs, (marginal_x * (ys[end] - ys[1]) / 25) .+ ys[1], legend=false)
                p = plot!(p, xs[end] .- (marginal_y * (xs[end] - xs[1]) / 25), ys)
                saveData && Export("$(PATH)$(title)_$(system)_$(i)_$(tstr)_marginal_x", xs, marginal_x)
                saveData && Export("$(PATH)$(title)_$(system)_$(i)_$(tstr)_marginal_p", ys, marginal_y)
            end

            if postProcess !== nothing
                p = postProcess(system, tout[i], p, postProcessParams)
            end

            if len > 0
                ps = plot(psp..., layout=operatorsLayout)
                p = plot(p, ps, layout=grid(1, 2, widths=[0.7,0.3]))
            end

            showGraph && display(plot(p))
            saveGraph && savefig(p, "$(PATH)$(title)_$(system)_$i$logstr.png")
        end

        println("display...$(time)s")

        saveData && Export("$(PATH)$(title)_$(system)_$tstr$logstr.txt", xs, ys, w)
    end

    return true
end

""" Complex survival amplitude """
function SurvivalAmplitude(system::QuantumSystem; Ψ0=nothing, mint=0.0, maxt=100.0, numt=1001, showGraph=true, saveGraph=true, saveData=true, kwargs...)    
    # Initial state
    if Ψ0 === nothing Ψ0 = ΨGS(system) end

    result = Array{Float64}(undef, 2, numt)

    time = @elapsed begin
        u = U(system, mint)
        du = U(system, (maxt - mint) / numt)
        Ψ = u * Ψ0      # Initial state

        for i in 1:numt
            s = sum(Ψ.data .* Ψ0.data)
            result[1, i] = real(s)
            result[2, i] = imag(s)
            Ψ = du * Ψ
        end
    end
    println(time)

    tout = LinRange(mint, maxt, numt)

    if showGraph || saveGraph
        ps = Array{Any}(undef, 4)
        ps[1] = plot(tout, result[1,:], xlabel="\$t\$", ylabel="Re\$A\$", legend=false)
        ps[2] = plot(tout, result[2,:], xlabel="\$t\$", ylabel="Im\$A\$", legend=false)
        ps[3] = plot(tout, log.(result[1,:].^2 .+ result[2,:].^2), xlabel="\$t\$", ylabel="\$\\ln|A|^2\$", legend=false)
        ps[4] = plot(result[1,:], result[2,:], xlabel="Re\$A\$", ylabel="Im\$A\$", legend=false)

        p = plot(ps..., layout=(2, 2), plot_title=String(system); kwargs...)
        showGraph && display(p)
        saveGraph && savefig(p, "$(PATH)amplitude_$system.pdf")
    end

    return result
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

        p = plot(ps..., layout=(2, (numop + 1) ÷ 2); kwargs...)
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

function ExpectationValuesLindblad(system::QuantumSystem, operators, lindblad; Ψ0=nothing, mint=0.0, maxt=100.0, numt=1001, showGraph=true, saveGraph=true, saveData=true, asymptotics=true, fname="", kwargs...)    
""" Expectation values of specific operator in Lindblad equation """
    # Initial state
    if Ψ0 === nothing Ψ0 = ΨGS(system) end

    println("Expectation values ", system, " E=", real(expect(H(system), Ψ0)) / Size(system), "...")

    numop = length(operators)
    result = Array{Float64}(undef, numop, numt)

    ts = LinRange(mint, maxt, numt)
    tout, ρt = timeevolution.master(ts, Ψ0, H(system), lindblad)

    time = @elapsed begin
        for (j, (name, operator)) in enumerate(operators)
            expectation = real(expect(operator, ρt))

            for i in 1:numt
                result[j, i] = expectation[i]
            end
        end
    end
    println("Finished ", system, " E=", real(expect(H(system), Ψ0)) / Size(system), " ", time)

    if asymptotics
        asymptoticValues = AsymptoticValues(system, operators; Ψ0=Ψ0)
    end

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
        saveGraph && savefig(p, "$(PATH)$(fname)lindblad_$system.png")
    end

    if saveData
        for (j, (name, _)) in enumerate(operators)
            Export("$(PATH)$(fname)$(name)_$system.txt", tout, result[j,:], asymptotics ? asymptoticValues[j] : nothing)
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
            print(".")
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