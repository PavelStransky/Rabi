using StatsBase
using DifferentialEquations
using LinearAlgebra

include("../../Rabi.jl")

const PATH = "d:/results/rabi/schnellbruder/"

gr()
default(size=(1000,1000))

function AlphaOnSphere(point; alpha_front=1.0, alpha_back=0.0)
    # Camera position
    camera_position = [5.0, -10.0, 5.0]

    # Vector from camera to point
    d = camera_position - point

    a = dot(d, d)
    b = 2 * dot(d, point)
    c = dot(point, point) - 1.0

    Δ = b^2 - 4*a*c
    if Δ < 0
        error("No real intersection with the sphere")
    end    

    t1 = (-b - sqrt(Δ)) / (2*a)
    t2 = (-b + sqrt(Δ)) / (2*a)

    t = max(t1, t2)

    return alpha_front - 0.5 * abs(t) * norm(d) * (alpha_front - alpha_back)
end


function PlotAxis(p)
    p = plot!(p, [0, 0], [0, 0], [1, 1.5], color=:black, lw=4)
    p = scatter!(p, [0], [0], [1], markersize=10, color=:black)
    p = scatter!(p, [0], [0], [-1], markersize=10, color=:black, alpha=AlphaOnSphere([0, 0, -1]))

    for z = -0.95:0.1:1
        p = plot!(p, [0, 0], [0, 0], [z, z + 0.04], color=:black, alpha=AlphaOnSphere([0, 0, z]), lw=4)
    end

    p = plot!(p, [0, 0], [0, 0], [-1.15, -1.05], color=:black, alpha=AlphaOnSphere([0, 0, -1.1]), lw=4)
    p = plot!(p, [0, 0], [0, 0], [-1.5, -1.15], color=:black, lw=4)
end


function PlotSphere()
    θ = range(0, stop=π, length=100)
    ϕ = range(0, stop=2π, length=100)
    x = [sin(t)*cos(p) for t in θ, p in ϕ]
    y = [sin(t)*sin(p) for t in θ, p in ϕ]
    z = [cos(t) for t in θ, p in ϕ]

    lim = 1.25

    p = surface(x, y, z, alpha=0.1, legend=false, color=:gray, camera=(30,30), 
    grid=false, ticks=nothing, axis=false,
    xlims=(-lim, lim), ylims=(-lim, lim), zlims=(-lim, lim))

    return p
end

" Determines whether the line is visible on the sphere or not"
function PlotLine(p, xs, ys, zs; color=:red, lw=2, alpha_front=0.8, alpha_back=0.0, marker=:circle)
    lastalpha = AlphaOnSphere([xs[1], ys[1], zs[1]])

    lx = []
    ly = []
    lz = []

    for (x, y, z) in zip(xs, ys, zs)
        alpha = AlphaOnSphere([x, y, z])

        push!(lx, x)
        push!(ly, y)
        push!(lz, z)

        if alpha != lastalpha
            p = plot!(p, lx, ly, lz, color=color, alpha=lastalpha, lw=lw)

            lastalpha = alpha
            lx = [x]
            ly = [y]
            lz = [z]
        end
    end

    p = plot!(p, lx, ly, lz, color=color, alpha=lastalpha, lw=lw)
    p = scatter!(p, [lx[end]], [ly[end]], [lz[end]], markersize=10, m=marker, color=color, alpha=lastalpha)

    return p
end

function Sphere(rabii::Rabi; λf=-0.37, min_time=0.0, max_time=300.0, num_time::Int=6000)	
    rabif = Copy(rabii, λ=λf)
    gs = SingleWellState(rabii)

    ev = ExpectationValues(rabif, [:Jx=>Jx(rabif), :Jy=>Jy(rabif), :Jz=>Jz(rabif)], Ψ0=gs, mint=min_time, maxt=max_time, numt=num_time, saveGraph=false, saveData=false, asymptotics=false)
    jx = ev[1, :] ./ rabif.j
    jy = ev[2, :] ./ rabif.j
    jz = ev[3, :] ./ rabif.j

    # n = sqrt.(jx .* jx .+ jy .* jy .+ jz .* jz)
    # jx = jx ./ n
    # jy = jy ./ n
    # jz = jz ./ n

    colours = palette(:auto)

    pa = PlotSphere()
          
    pa = PlotLine(pa, jx, jy, jz, color=:red, lw=2, marker=:diamond)

    pa = PlotAxis(pa)
    display(pa)

    savefig(pa, PATH * "sphere_$(min_time).png")
end

Sphere(Rabi(R=50, λ=1.5, δ=0.5, j=1//2); λf=-24/65, min_time=0.0, max_time=300.0, num_time=6000)