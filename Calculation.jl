using Distributed

if nprocs() <= WORKERS
    addprocs(WORKERS + 1 - nprocs())
end

@everywhere include("Export.jl")
@everywhere include("Rabi.jl")
@everywhere include("QuarticOscillator.jl")

@everywhere ComputeAsymptotics(rabi::Rabi) = AsymptoticValuesMatrix(rabi, AllOperators(rabi))
@everywhere ComputeTime(rabi::Rabi; mint=0, maxt=20, numt=201) = ExpectationValues(rabi, AllOperators(rabi); mint=mint, maxt=maxt, numt=numt, showGraph=false, saveGraph=false, saveData=false)

function Plot(xs, result)
    p1 = plot(xs, result[:,1], title="Jx", ylims=(-0.5, 0.5))
    p2 = plot(xs, result[:,2], title="Jy", ylims=(-0.5, 0.5))
    p3 = plot(xs, result[:,3], title="Jz", ylims=(-0.5, 0.5))

    p4 = plot(xs, result[:,4], title="P", ylims=(0, 1))
    p5 = plot(xs, result[:,5], title="q")
    p7 = plot(xs, result[:,7], title="n")

    p = plot(p1, p2, p3, p4, p5, p7, layout=(3, 2))
    display(p)    
end

function Compute(rabi::Rabi; λs=nothing, δs=nothing, Rs=nothing, parallel=true)  
    isλ = λs !== nothing
    isδ = δs !== nothing
    isR = Rs !== nothing

    result = nothing
    xs = nothing

    startTime = time()
    if isλ && isδ
        input = [[λ, δ] for λ in λs, δ in δs]
        if parallel
            result = pmap((args)->ComputeAsymptotics(Rabi(rabi=rabi, N=-1, λ=args[1], δ=args[2])), input)
            result = [result[i,j][k] for k in 1:length(result[1,1]), i in 1:length(λs), j in 1:length(δs)]
        end
    else
        if isλ
            result = parallel ? pmap((λ)->ComputeAsymptotics(Rabi(rabi; N=-1, λ=λ)), λs) : [ComputeAsymptotics(Rabi(rabi; N=-1, λ=λ)) for λ in λs]
            xs = λs
        elseif isδ
            result = parallel ? pmap((δ)->ComputeAsymptotics(Rabi(rabi; N=-1, δ=δ)), δs) : [ComputeAsymptotics(Rabi(rabi; N=-1, δ=δ)) for δ in δs]
            xs = δs
        elseif isR
            result = parallel ? pmap((R)->ComputeAsymptotics(Rabi(rabi; N=-1, R=R)), Rs) : [ComputeAsymptotics(Rabi(rabi; N=-1, R=R)) for R in Rs]
            xs = Rs
        end
        result = [result[i][j] for i in 1:length(result), j in 1:length(result[1])]
        Plot(xs, result)
    end
    println(time() - startTime)

    fname = string(PATH, "map",
        isλ ? "L" : "", 
        isδ ? "D" : "", 
        isR ? "R" : "",
        "_(",
        isR ? "" : "$(rabi.R),", 
        isλ ? "" : "$(rabi.λ),", 
        isδ ? "" : "$(rabi.δ),",
        "$(rabi.μ),$(rabi.ν))")
    
    if isλ && isδ    
        Export(fname, λs, δs, result)
    else
        Export(fname, xs, result)
    end

    return result
end

function DQPT(rabi::Rabi; λs=nothing, μs=nothing, mint=0.0, maxt=30.0, numt=601, showGraph=true, parallel=true)
    isλ = λs !== nothing
    isμ = μs !== nothing

    time = @elapsed begin
        if isλ
            if parallel
                result = pmap((λ)->ComputeTime(Copy(rabi, N=-1, λ=λ); mint=mint, maxt=maxt, numt=numt), λs)
            else
                result = [ComputeTime(Copy(rabi, N=-1, λ=λ); mint=mint, maxt=maxt, numt=numt) for λ in λs]
            end
            xs = λs
            fname = "dqptL_($(rabi.R),$(rabi.δ),$(rabi.μ),$(rabi.ν))"

        elseif isμ
            if parallel
                result = pmap((μ)->ComputeTime(Copy(rabi, N=-1, μ=μ, ν=μ); mint=mint, maxt=maxt, numt=numt), μs)
            else
                result = [ComputeTime(Copy(rabi, N=-1, μ=μ, ν=μ); mint=mint, maxt=maxt, numt=numt) for μ in μs]
            end
            
            xs = log10.(μs)
            fname = "dqptM_($(rabi.R),$(rabi.λ),$(rabi.δ))"
        end

        # Reshape of the results
        result = [result[k][i,j] for i in 1:length(result[1][:,1]), k in 1:length(xs), j in 1:numt]

    end
    println(time)

    ts = LinRange(mint, maxt, numt)
    Export("$(PATH)$fname", ts, xs, result)

    if showGraph
        p1 = heatmap(ts, xs, result[1,:,:], color=:coolwarm, clims=(-0.1, 0.1))
        p2 = heatmap(ts, xs, result[2,:,:], color=:coolwarm, clims=(-0.1, 0.1))
        p3 = heatmap(ts, xs, result[3,:,:], color=:coolwarm, clims=(-0.5, 0.5))
        p4 = heatmap(ts, xs, result[4,:,:], clims=(0, 1))
        p5 = heatmap(ts, xs, result[5,:,:], color=:coolwarm, clims=(-0.5, 0.5))
        p6 = heatmap(ts, xs, result[6,:,:], color=:coolwarm, clims=(-0.5, 0.5))
    
        display(plot(p1, p2, p3, p4, p5, p6))
    end

    return result
end
