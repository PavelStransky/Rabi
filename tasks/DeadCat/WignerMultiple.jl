using Distributed

@everywhere WORKERS = 1
include("../../Calculation.jl")
@everywhere const PATH = "d:/results/Rabi/deadcat/multiple/"
@everywhere pyplot(size = (1200, 1200))

function WignerMultiple(systems, tss; xs=LinRange(-1, 1, 101), ys=nothing, clim=(-0.2, 0.2), title="", kwargs...)
    """ Wigner function for four different Rabi models"""

    # Range in y direction
    if ys === nothing ys = xs end
    
    numt = length(tss[1])

    for i = 67:numt
        ps = Array{Any}(undef, length(systems))

        for (index, system) = enumerate(systems)
            # Initial state
            Ψ0 = ΨGS(system)

            ts = [0, tss[index][i]]

            # Schroedinger time evolution
            print("Schrödinger ", system, "...")
            time = @elapsed tout, Ψt = timeevolution.schroedinger(ts, Ψ0, H(system))
            println(time, "s")

            xss = xs * sqrt(2.0 * Size(system))
            yss = ys * sqrt(2.0 * Size(system))

            tstr = format(tout[end], precision=2)

            print("T[$i]=$tstr trace...")
            time = @elapsed ρ = DensityMatrix(system, Ψt[2])
            print("$time $title...")

            time = @elapsed w = wigner(ρ, xss, yss)
            print("$time display...")

            time = @elapsed begin
                p = heatmap(xs, ys, transpose(w), color=:bwr, grid=false, title="\$t=$tstr\$", xlabel=raw"$q$", ylabel=raw"$p$", legend=false, clim=clim, kwargs...)

                marginal_x = vec(sum(w, dims=2))
                marginal_x = marginal_x / (sum(marginal_x) * (xs[2] - xs[1]))
                marginal_y = vec(sum(w, dims=1))
                marginal_y = marginal_y / (sum(marginal_y) * (ys[2] - ys[1]))
                p = plot!(p, xs, (marginal_x * (ys[end] - ys[1]) / 15) .+ ys[1], legend=false)
                p = plot!(p, xs[end] .- (marginal_y * (xs[end] - xs[1]) / 15), ys)
                ps[index] = p
            end
        end

        p = plot(ps..., layout=grid(2, 2))
        # display(p)
        savefig(p, "$(PATH)$(title)_$i.png")    

    end
    
end

function Run()
    λ = 0.75
    δ = 0.5

    # Final values for the figure
    t = LinRange(0, 2.5, 126)
    ts = 20 * t

    R = 25
    μ = 0.132 / R
    ν = μ
    rabi1 = Rabi(R=R, δ=δ, λ=λ, μ=μ, ν=ν)
    t1 = 17.8 * t 
    ts1 = ts

    R = 50
    μ = 0.132 / R
    ν = μ
    rabi2 = Rabi(R=R, δ=δ, λ=λ, μ=μ, ν=ν)
    t2 = 18.88 * t
    ts2 = 1.05564755 .+ ts

    R = 100
    μ = 0.132 / R
    ν = μ
    rabi3 = Rabi(R=R, δ=δ, λ=λ, μ=μ, ν=ν)
    t3 = 20.31 * t
    ts3 = 2.48564755 .+ ts

    R = 200
    μ = 0.132 / R
    ν = μ
    rabi4 = Rabi(R=R, δ=δ, λ=λ, μ=μ, ν=ν)
    t4 = 21.71 * t
    ts4 = 3.88564755 .+ ts

    num = 401

     wigner = WignerMultiple([rabi1, rabi2, rabi3, rabi4], [ts, ts, ts, ts]; xs=LinRange(-1.2, 1.2, num), ys=LinRange(-0.6, 0.6, num), clim=(-0.1, 0.1), title="linear")
#    wigner = WignerMultiple([rabi1, rabi2, rabi3, rabi4], [t1, t2, t3, t4]; xs=LinRange(-1.2, 1.2, num), ys=LinRange(-0.6, 0.6, num), clim=(-0.1, 0.1), title="scaled")
#    wigner = WignerMultiple([rabi1, rabi2, rabi3, rabi4], [ts1, ts2, ts3, ts4]; xs=LinRange(-1.2, 1.2, num), ys=LinRange(-0.6, 0.6, num), clim=(-0.1, 0.1), title="shifted")
end

Run()
exit()