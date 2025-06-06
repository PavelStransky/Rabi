using StatsBase
using DifferentialEquations
using LinearAlgebra

# Run at IPNP36:
#     - python3 CreateJobs.py
#       creates jobs.txt file with all the jobs
#     - parallel -j 24 --eta --joblog joblog.txt --results logs < jobs.txt
#       runs jobs.txt at 24 cores using a simple "parallel" scheduler; jobs are "nice" and all results are logged

#     Creating the video:
#     ffmpeg -f lavfi -t 3 -i color=black:s=1920x1080:r=30 -framerate 30 -i "%d.png" -i music.mp3 -filter_complex "[1:v]scale=iw*min(1920/iw\,1080/ih):ih*min(1920/iw\,1080/ih),pad=1920:1080:(1920-iw*min(1920/iw\,1080/ih))/2:(1080-ih*min(1920/iw\,1080/ih))/2:color=black,fade=t=in:st=0:d=1[img]; [0:v][img]concat=n=2:v=1:a=0[outv]" -map "[outv]" -map 2:a -preset slow -crf 18 -c:v libx264 -pix_fmt yuv420p -c:a aac -profile:v high out.mp4

#     With title (lf=0.5):
#     ffmpeg -framerate 30 -loop 1 -t 9 -i intro.png -framerate 30 -i "%d.png" -i music.mp3 -filter_complex "[1:v]scale=iw*min(1920/iw\,1080/ih):ih*min(1920/iw\,1080/ih),pad=1920:1080:(1920-iw*min(1920/iw\,1080/ih))/2:(1080-ih*min(1920/iw\,1080/ih))/2:color=black,setsar=1,fade=t=in:st=0:d=1[img]; [0:v]setsar=1[title]; [title][img]concat=n=2:v=1:a=0[outv]" -map "[outv]" -map 2:a -preset slow -crf 18 -c:v libx264 -pix_fmt yuv420p -c:a aac -profile:v high out.mp4

#     Without audio:
#     ffmpeg -f lavfi -t 1 -i color=black:s=1920x1080:r=30 -framerate 30 -i "%d.png" -filter_complex "[1:v]scale=iw*min(1920/iw\,1080/ih):ih*min(1920/iw\,1080/ih),pad=1920:1080:(1920-iw*min(1920/iw\,1080/ih))/2:(1080-ih*min(1920/iw\,1080/ih))/2:color=black,fade=t=in:st=0:d=1[img]; [0:v][img]concat=n=2:v=1:a=0[outv]" -map "[outv]" -preset slow -crf 18 -c:v libx264 -pix_fmt yuv420p -profile:v high out.mp4

#Suppress strange warnings if using pyplot
# using PyCall
# pyimport("warnings").filterwarnings("ignore", message=".*No data for colormapping provided.*")

include("../Rabi.jl")

const PATH = "d:/results/rabi/schnellbruder/"
# const PATH = "/home/stransky/results/"

gr()
default(size=(1920,1080))

""" Evolution of the partial trace and related quantities """
function PartialTraceEvolution(rabii; λf=0.5, mint=0.0, maxt=50, numt=1000, log=false)
    # Initial and final systems
    rabif = Copy(rabii; λ=λf)
    gs = SingleWellState(rabii)

    spin = Array{Float64}(undef, numt, 4)
    ptt = Array{Float64}(undef, numt, 4)
    ptev = Array{Float64}(undef, numt, Int(2 * j + 1))
    ptdiag = Array{Float64}(undef, numt, Int(2 * j + 1))
    field = Array{Float64}(undef, numt, 4)

    time = @elapsed begin
        u = U(rabif, mint)
        du = U(rabif, (maxt - mint) / numt)
        Ψ = u * gs      # Initial state

        for i in 1:numt
            ρq = ptrace(Ψ, 1)
            es, vs = QuantumOptics.eigenstates(ρq)

            for j in 1:Int(2 * j + 1)
                ptev[i, j] = if log log10(real(es[j])) else real(es[j]) end
                ptdiag[i, 1] = real(ρq.data[j, j])
            end

            ptt[i, 1] = real(ρq.data[1, 2])
            ptt[i, 2] = -imag(ρq.data[1, 2])
            ptt[i, 3] = real(ρq.data[1, 1]) - 0.5
            ptt[i, 4] = sqrt(ptt[i, 1]^2 + ptt[i, 2]^2 + ptt[i, 3]^2)

            q = ExpectationValue("X", X(rabif), Ψ, rabif)
            p = ExpectationValue("P", P(rabif), Ψ, rabif)

            field[i, 1] = sqrt(8) * rabif.λ * q
            field[i, 2] = -sqrt(8) * rabif.λ * rabif.δ * p
            field[i, 3] = 1 / (2 * rabif.j) + sqrt(8) * rabif.μ * q
            field[i, 4] = sqrt(field[i, 1]^2 + field[i, 2]^2 + field[i, 3]^2)

            for j in 1:3
                field[i, j] = -field[i, j] / field[i, 4]
            end

            spin[i, 1] = ExpectationValue("Jx", Jx(rabif), Ψ, rabif)
            spin[i, 2] = ExpectationValue("Jy", Jy(rabif), Ψ, rabif)
            spin[i, 3] = ExpectationValue("Jz", Jz(rabif), Ψ, rabif)
            spin[i, 4] = sqrt(spin[i, 1]^2 + spin[i, 2]^2 + spin[i, 3]^2)

            for j in 1:3
                spin[i, j] = spin[i, j] / spin[i, 4]
            end

            Ψ = du * Ψ
        end
    end
    println(time)

    tout = LinRange(mint, maxt, numt)

    p = plot(tout, field[:, 1:3], xlabel="t", lc = [:red :blue :green], label=["bx" "by" "bz"], lw=2, alpha=0.5)
    p = plot!(p, tout, spin[:, 1:3], lc = [:red :blue :green], label=["jx" "jy" "jz"], lw=1, ls=:dash, alpha=0.5)
    display(p)
    savefig(p, "$(PATH)components_($(rabif)).png")

    if log
        p = plot(tout, ptev, ylims=(-5, 0), lc = :black, xlabel="\$t\$", label=["ρ_spin EV 1" "ρ_spin EV 2"], ls=[:dash :solid])
        q = plot(tout, ptdiag, ylims=(-5, 0), lc = :black, xlabel="\$t\$", label=["ρ_diag 1" "ρ_diag 2"], ls=[:dash :solid])
    else
        p = plot(tout, ptev, lc = :black, xlabel="\$t\$", label=["ρ_spin EV 1" "ρ_spin EV 2"], ls=[:dash :solid])
        q = plot(tout, ptdiag, lc = :black, xlabel="\$t\$", label=["ρ_diag 1" "ρ_diag 2"], ls=[:dash :solid])
    end

    display(p)
    savefig(p, "$(PATH)ptt_($(rabif)).png")
end

""" Strength function for the Rabi model and its local extrema """
function StrengthFunction(rabii; λf=0.5, threshold=1E-15, envelope_window=3, showgraph=true, savedata=false)
    """ Finds all local extrema in the array a """
    function LocalExtremaIndices(a)
        n = length(a)
        extrema_indices = Int[]

        # Include the first...
        push!(extrema_indices, 1)

        for i in 2:(n-1)
            if (a[i] > a[i-1] && a[i] > a[i+1]) || (a[i] < a[i-1] && a[i] < a[i+1])
                push!(extrema_indices, i)
            end
        end
    
        # ...and the last elements as extrema
        push!(extrema_indices, n)

        return extrema_indices
    end    

    # Initial and final systems
    rabif = Copy(rabii; λ=λf)
    gs = SingleWellState(rabii)

    # Strength function
    sf = StrengthFunction(rabif, gs)

    energies = sf[1]
    probabilities = sf[2]

    indices = findall(p -> p > threshold, probabilities)

    energies = energies[indices]
    probabilities = probabilities[indices]

    envelopex = copy(energies)
    envelopey = copy(probabilities)

    for i in 1:length(probabilities)
        k = argmax(probabilities[max(i - envelope_window, 1):i]) + max(i - envelope_window, 1) - 1
        envelopex[i] = energies[k]
        envelopey[i] = probabilities[k]
    end

        # Find indices of the first occurrences of unique y-values
    uniquey = unique(envelopey)

    unique_indices = Int[]
    for y in uniquey
        push!(unique_indices, findfirst(e -> e == y, envelopey))
    end

    # Filter x and y arrays based on these indices
    envelopex = envelopex[unique_indices]
    envelopey = envelopey[unique_indices]

    extrema_indices = LocalExtremaIndices(envelopey)
    
    println("Local extrema indices: ", extrema_indices)
    println("Local extrema x: ", envelopex[extrema_indices])
    println("Local extrema y: ", envelopey[extrema_indices])

    sum_peaks = Array{Float64}(undef, Int(2 * rabii.j + 1))

    k = 1
    for i in 3:length(extrema_indices)
        if isodd(i)
            sum_peaks[k] = sum(envelopey[findall(e -> envelopex[extrema_indices[i - 2]] < e < envelopex[extrema_indices[i]], envelopex)])
            println("Peak $(k): ", sum_peaks[k])
            k += 1
        end
    end

    println("Peak check: ", sum(sum_peaks))

    if showgraph
        p = scatter(sf[1], log10.(sf[2]),  ylims=(log10(threshold) - 5, 0), xlims=(0,3), title="$(k)", markeralpha=0.5, markerstrokewidth=0, xlabel="\$E\$", ylabel=raw"$\log_{10}S$", label=raw"$|\psi_{\mathrm{GS}}\rangle$")
        p = scatter!(p, envelopex, log10.(envelopey), color=:red, markerstrokewidth=0)
        display(plot(p))
    end
    if savedata
        Export("$(PATH)sf_$(String(rabii))_$(λf)", sf[1], sf[2])
        Export("$(PATH)sf_$(String(rabii))_$(λf)_envelope", envelopex, envelopey)
    end

    return sum_peaks
end

function EquationOfMotion!(dx, x, parameters, t)
    q, p = x
    rabi, m = parameters

    s2 = 2 * rabi.λ^2 * (q^2 + rabi.δ^2 * p^2) + 1    
    s = sqrt(s2)

    dx[1] = p * (1 + rabi.λ^2 * rabi.δ^2 * m / rabi.j / s)
    dx[2] = -q * (1 + rabi.λ^2 * m / rabi.j / s)
end

function ClassicalToWigner(rabif, t, p, λi)
    x0 = [-sqrt(0.5 * (λi^2 - 1 / λi^2)), 0]
    
    timeInterval = (0.0, t)
    solver = TsitPap8()
    tolerance = 1E-8
    fnc = ODEFunction(EquationOfMotion!)

    colours = palette(:auto)
    colourIndex = 3

    for m = -rabif.j:rabif.j
        problem = ODEProblem(fnc, x0, timeInterval, (rabif, m))

        saveat = collect(range(max(t - 10, 0), t, step=0.05))
        time = @elapsed solution = solve(problem, solver, reltol=tolerance, abstol=tolerance, verbose=true, saveat=saveat)

        println("m = $m, time = $(time)s, solution = $(length(solution))")

        p = plot!(p, solution, idxs=(1,2), lw=2, alpha=0.4, c=colours[colourIndex], label="m=$m", xlabel=raw"$q$", ylabel=raw"$p$")
        p = scatter!(p, [solution(t)[1]], [solution(t)[2]], markersize=10, c=colours[colourIndex], label=nothing)

        colourIndex += 1
    end

    return p
end

function WignerFunctions(rabii; λf=0.5, wignerMesh=301, limits=1.6, maxt=30, numt=31, log=false, firstIndex=1, lastIndex=-1, kwargs...)
    # Initial and final systems
    rabif = Copy(rabii; λ=λf)
    gs = SingleWellState(rabii)	

    if log
        clim = (-6, 6)
    else
        clim = (-0.08, 0.08)
    end

    Wigner(rabif; Ψ0=gs, operators=[:Jx=>Jx(rabif), :q=>X(rabif), :Jy=>Jy(rabif), :p=>P(rabif), :Jz=>Jz(rabif), :P=>projector(gs)], operatorsColor=[:gray, :red, :gray, :red, :gray, :blue], operatorsMarker=[:diamond, :square, :diamond, :square, :diamond, :pentagon],
    operatorsLayout=(7, 2), firstIndex=firstIndex, lastIndex=lastIndex,
    maxt=maxt, numt=numt, xs=LinRange(-limits, limits, wignerMesh), ys=LinRange(-limits, limits, wignerMesh), 
    clim=clim, saveData=false, saveGraph=true, showGraph=true, log=log, postProcess=ClassicalToWigner, postProcessParams=rabii.λ, kwargs...)
end

function PartialTraceEvolutionScaling(rabii; R0=10, numR=6, λi=1.5, mint=0.0, maxt=50, numt=1000)
    pt = Array{Float64}(undef, numt, 2 * numR)

    tout = LinRange(mint, maxt, numt)
    lc = []
    label = []

    colors = [:red :blue :green :purple :orange :black]

    R = R0
    for j = 1:numR
        print(R)
        print("...")

        # Initial and final systems
        rabif = Copy(rabii; λ=λf)
        gs = SingleWellState(rabii)
    
        time = @elapsed begin
            u = U(rabif, mint)
            du = U(rabif, (maxt - mint) / numt)
            Ψ = u * gs      # Initial state

            for i in 1:numt
                ρq = ptrace(Ψ, 1)
                es, vs = QuantumOptics.eigenstates(ρq)

                pt[i, 2*j - 1] = real(es[1])
                pt[i, 2*j] = real(es[2])

                Ψ = du * Ψ
            end
        end
        println(time)

        push!(lc, colors[j])
        push!(lc, colors[j])

        push!(label, "EV1 R=$R")
        push!(label, "EV2 R=$R")

        R *= 2
    end

    lc = permutedims(lc)
    label = permutedims(label)

    p = plot(tout, pt, xlabel="\$t\$", lc=lc, label=label)
    display(p)
    savefig(p, "$(PATH)ptt_scaling.png")
end

# PartialTraceEvolution(λf=0.0, maxt=1.0)
# WignerFunctions(λf=0.0)

# PartialTraceEvolution(λf=0.5)
# WignerFunctions(λf=0.5)

# PartialTraceEvolution(λf=1.0, maxt=100)
# WignerFunctions(λf=1.0, log=true)

# PartialTraceEvolutionScaling()
# PartialTraceEvolutionScaling(λf=1.0, maxt=200)

# PartialTraceEvolution(λf=0.5, j=4//2, log=true, maxt=100)
# WignerFunctions(λf=0.5, j=4//2, log=true, range=1.8, ts=LinRange(0, 100, 51), wignerMesh=401)

# PartialTraceEvolution(λf=0.5, δ=0.0)

function alpha2(λi, λf)
    return ((λf - 16 * λf * λi^4 + λi * (-1 + 4 * sqrt(λi^4 + λf^2 * λi^2 * (-1 + 16 * λi^4)))) /
              (λi * sqrt(λi^4 + λf^2 * λi^2 * (-1 + 16 * λi^4)))) / 8
end

function beta2(λi, λf)
    return 0.25 * (2 + (λi + λf * (-1 + 16 * λi^4)) / (2 * λi * sqrt(λi^4 + λf^2 * λi^2 * (-1 + 16 * λi^4))))
end

function StrengthFunctionMagnitude(rabii; λf_min=-0.5, λf_max=1.0, λf_num=50)
    peaks = Array{Float64}(undef, λf_num, Int(2 * rabii.j + 1))
    λfs = LinRange(λf_min, λf_max, λf_num)

    for i = eachindex(λfs)
        λf = λfs[i]
        ps = StrengthFunction(rabii; λf=λf, showgraph=false, envelope_window=Int(6 * rabii.j + 2))

        for k = eachindex(ps)
            peaks[i, k] = ps[k]
        end
    end

    p = plot(λfs, peaks, xlabel="\$λ_f\$", ylabel="Strength", title="Strength function")
    display(p)

    return λfs, peaks
end

# _, pld = LevelDynamics(Rabi(R=10,N=300), ps=LinRange(-1.5, 1.5, 401), limit=150, ylims=(-1,2), saveGraph=false)

if length(ARGS) > 0
    firstIndex = 30 * parse(Int, ARGS[1]) + 1
    lastIndex = firstIndex + 29
    println("firstIndex = $(firstIndex), lastIndex = $(lastIndex)")
    println(Plots.backend())

    type = 0
    if length(ARGS) > 1
        type = parse(Int, ARGS[2])
    end

    λf = -0.37
    rabi = Rabi(R=50, λ=1.5, δ=0.5, j=2//2)

    maxt = 300
    numt = 6000
    limits = 1.5

    if type == 1
        rabi = Rabi(R=50, λ=1.5, δ=0.0, j=2//2)
    elseif type == 2
        rabi = Rabi(R=50, λ=1.5, δ=0.5, j=4//2)
    elseif type == 3
        rabi = Rabi(R=50, λ=1.5, δ=0.0, j=4//2)
    elseif type == 4
        λf = -1.5
        maxt = 120
        numt = 6000
        rabi = Rabi(R=50, λ=1.5, δ=0.5, j=1//2)
    elseif type == 5
        λf = -1.0
        maxt = 120
        numt = 6000
        rabi = Rabi(R=50, λ=1.5, δ=0.5, j=1//2)
    elseif type == 6
        λf = 0.5
        maxt = 120
        numt = 6000
        rabi = Rabi(R=50, λ=1.5, δ=0.5, j=1//2)
    elseif type == 7
        λf = -sqrt(209) / 12
        maxt = 120
        numt = 6000
        rabi = Rabi(R=50, λ=1.5, δ=0.5, j=1//2)
    elseif type == 8
        maxt = 120
        numt = 6000
        rabi = Rabi(R=50, λ=1.5, δ=0.5, j=1//2)
    elseif type == 9
        λf = -2.0
        limits = 2.0
        maxt = 60
        numt = 6000
        rabi = Rabi(R=50, λ=1.5, δ=0.5, j=1//2)
    elseif type == 10
        λf = -1.5
        rabi = Rabi(R=50, λ=1.5, δ=0.5, j=1//2)
    elseif type == 11
        λf = -1.0
        rabi = Rabi(R=50, λ=1.5, δ=0.5, j=1//2)
    elseif type == 12
        λf = 0.5
        rabi = Rabi(R=50, λ=1.5, δ=0.5, j=1//2)
    elseif type == 13
        λf = -1.5
        rabi = Rabi(R=50, λ=1.5, j=1//2)
    elseif type == 14
        λf = -1.0
        rabi = Rabi(R=50, λ=1.5, j=1//2)
    elseif type == 15
        λf = 0.5
        rabi = Rabi(R=50, λ=1.5, j=1//2)
    elseif type == 20
        λf = -0.62116
        maxt = 120
        numt = 6000
        rabi = Rabi(R=50, λ=1.5, δ=0.5, j=1//2)
    elseif type == 21
        λf = 0.62116
        maxt = 120
        numt = 6000
        rabi = Rabi(R=50, λ=1.5, δ=0.5, j=1//2)
    end

    WignerFunctions(rabi, λf=λf, limits=limits, wignerMesh=501, maxt=maxt, numt=numt, showGraph=false, firstIndex=firstIndex, lastIndex=lastIndex, marginals=true)

    exit()
end

function Figure1()
    default(size=(800,600))

    rabi = Rabi(R=50, λ=1.5, δ=0.5, j=1//2)
    gs = SingleWellState(rabi)

    ef, sf, parity = StrengthFunction(Copy(rabi; λ=-0.3), gs)

    # "Odd" and "even" parity states
    sfp = Array{Float64}(undef, 0)
    ep = Array{Float64}(undef, 0)
    sfn = Array{Float64}(undef, 0)
    en = Array{Float64}(undef, 0)
    for (j, p) in enumerate(parity)
        if p > 0
            append!(ep, ef[j])
            append!(sfp, sf[j])
        else
            append!(en, ef[j])
            append!(sfn, sf[j])
        end
    end

    p = scatter(ep, sfp, markeralpha=0.9, markerstrokewidth=0, markersize=5, xlims=(0,3), label="Parity +", xlabel="\$e\$", ylabel="\$S\$", title=String(rabi))
    p = scatter(p, en, sfn, markeralpha=0.9, markerstrokewidth=0, markersize=5, label="Parity -")
    # p = plot(p, [0, λ], [-0.5, ExpectationValue(H(rabi), Ψ) / rabi.R], linewidth=3, arrow=(:closed, 1.0))
    display(plot(p))
    savefig(p, "d:/sf.png")
    # print(pf)

    # p = scatter(ef, sf, markeralpha=0.5, markerstrokewidth=0, xlabel="\$E\$", ylabel=raw"$\log_{10}S$", label=raw"$|\psi_{\mathrm{GS}}\rangle$")
    # display(plot(p))
end

# rabi = Rabi(R=50, λ=1.5, δ=0.5, j=1//2)
# x, y = StrengthFunctionMagnitude(rabi; λf_min=-1.2, λf_max=1.2, λf_num=200)
# Export("$(PATH)sf_$(String(rabi))", x, y)

# rabi = Rabi(R=50, λ=1.5, δ=0.5, j=2//2)
# x, y = StrengthFunctionMagnitude(rabi; λf_min=-1.2, λf_max=1.2, λf_num=200)
# Export("$(PATH)sf_$(String(rabi))", x, y)

# rabi = Rabi(R=50, λ=1.5, δ=0.5, j=4//2)
# x, y = StrengthFunctionMagnitude(rabi; λf_min=-1.2, λf_max=1.2, λf_num=200)
# Export("$(PATH)sf_$(String(rabi))", x, y)

# ps = StrengthFunction(rabi; λf=-1.0, showgraph=true, envelope_window=Int(6 * rabi.j + 2))


# Strength function
# Figure1()
rabi = Rabi(R=50, λ=1.5, δ=0.5, j=1//2)
StrengthFunction(rabi; λf=1.20474, showgraph=true, savedata=true, envelope_window=Int(6 * rabi.j + 2), threshold=1E-10)
StrengthFunction(rabi; λf=0.5, showgraph=true, savedata=true, envelope_window=Int(6 * rabi.j + 2), threshold=1E-10)
StrengthFunction(rabi; λf=0.0, showgraph=true, savedata=true, envelope_window=Int(6 * rabi.j + 2), threshold=1E-10)
StrengthFunction(rabi; λf=-0.369, showgraph=true, savedata=true, envelope_window=Int(6 * rabi.j + 2), threshold=1E-10)
StrengthFunction(rabi; λf=-1.20474, showgraph=true, savedata=true, envelope_window=Int(6 * rabi.j + 2), threshold=1E-10)

rabi = Rabi(R=50, λ=1.5, δ=0.5, j=4//2)
StrengthFunction(rabi; λf=-0.92678, showgraph=true, savedata=true, envelope_window=Int(6 * rabi.j + 2), threshold=1E-10)
StrengthFunction(rabi; λf=-0.57987, showgraph=true, savedata=true, envelope_window=Int(6 * rabi.j + 2), threshold=1E-10)
StrengthFunction(rabi; λf=-0.369, showgraph=true, savedata=true, envelope_window=Int(6 * rabi.j + 2), threshold=1E-10)
StrengthFunction(rabi; λf=-0.19734, showgraph=true, savedata=true, envelope_window=Int(6 * rabi.j + 2), threshold=1E-10)
StrengthFunction(rabi; λf=-0.02001, showgraph=true, savedata=true, envelope_window=Int(6 * rabi.j + 2), threshold=1E-10)
