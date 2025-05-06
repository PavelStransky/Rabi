using StatsBase
using DifferentialEquations
using LinearAlgebra

include("../Rabi.jl")

const PATH = "d:/results/rabi/schnellbruder/"
# const PATH = "/home/stransky/results/"

gr()
default(size=(1000,1000))


function EquationOfMotion!(dx, x, parameters, t)
    q, p = x
    rabi, m = parameters

    s2 = 2 * rabi.λ^2 * (q^2 + rabi.δ^2 * p^2) + 1   
    s = sqrt(s2)

    dx[1] = p * (1 + rabi.λ^2 * rabi.δ^2 * m / rabi.j / s)
    dx[2] = -q * (1 + rabi.λ^2 * m / rabi.j / s)
end

function Trajectory(rabif, λi, max_time, num_time)
    # Initial state
    x0 = [-sqrt(0.5 * (λi^2 - 1 / λi^2)), 0]

    timeInterval = (0.0, max_time)
    solver = TsitPap8()
    tolerance = 1E-8
    fnc = ODEFunction(EquationOfMotion!)

    result = []

    for m = -rabif.j:rabif.j
        problem = ODEProblem(fnc, x0, timeInterval, (rabif, m))

        saveat = collect(range(0, max_time, step=max_time / num_time))
        solution = solve(problem, solver, reltol=tolerance, abstol=tolerance, verbose=true, saveat=saveat)

        if length(result) == 0
            push!(result, solution.t)
        end

        push!(result, [v[1] for v in solution.u])
        push!(result, [v[2] for v in solution.u])
    end

    return result
end


function VisibilityOnSphere(point)
    if norm(point) < 0.99
        return false
    end

    # Camera position
    camera_position = [5.0, -10.0, 5.0]

    # Vector from camera to point
    camera_to_point = point - camera_position

    # Find the distance from the camera to the point along the line of sight
    d_camera_to_sphere = norm(camera_position)

    # Check if the point is hidden behind the sphere
    # If the point's distance to the camera is greater than the distance to the sphere's surface
    # then the point is visible (in front); otherwise, it’s hidden (behind).
    distance_to_point = norm(camera_to_point)

    # Determine visibility
    return distance_to_point < d_camera_to_sphere
end


function PlotAxis(p)
    p = plot!(p, [0, 0], [0, 0], [1, 1.5], color=:black, lw=4)
    p = scatter!(p, [0], [0], [1], markersize=10, color=:black)
    p = scatter!(p, [0], [0], [-1], markersize=10, color=:black, alpha=0.3)

    for z = -0.95:0.1:1
        p = plot!(p, [0, 0], [0, 0], [z, z + 0.04], color=:black, alpha=0.3, lw=4)
    end

    p = plot!(p, [0, 0], [0, 0], [-1.15, -1.05], color=:black, alpha=0.3, lw=4)
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
function PlotLine(p, xs, ys, zs; color=:red, lw=2, alpha_front=0.7, alpha_back=0.2, marker=:circle)
    visible = VisibilityOnSphere([xs[1], ys[1], zs[1]])

    lx = []
    ly = []
    lz = []

    for (x, y, z) in zip(xs, ys, zs)
        v = VisibilityOnSphere([x, y, z])

        push!(lx, x)
        push!(ly, y)
        push!(lz, z)

        if v != visible
            p = plot!(p, lx, ly, lz, color=color, alpha=visible ? alpha_front : alpha_back, lw=lw)
            
            visible = v
            lx = [x]
            ly = [y]
            lz = [z]
        end
    end

    p = plot!(p, lx, ly, lz, color=color, alpha=visible ? alpha_front : alpha_back, lw=lw)
    p = scatter!(p, [lx[end]], [ly[end]], [lz[end]], markersize=10, m=marker, color=color, alpha=visible ? 1.0 : 0.5)

    return p
end

function SphereAnimation(rabii::Rabi; λf=-0.37, max_time=300, num_time::Int=6000, show_time::Int=5, ev_coef=20)	
    rabif = Copy(rabii, λ=λf)
    gs = SingleWellState(rabii)

    trajectory = Trajectory(rabif, rabii.λ, max_time, num_time)

    ev = ExpectationValues(rabif, [:Jx=>Jx(rabif), :Jy=>Jy(rabif), :Jz=>Jz(rabif), :p=>P(rabif), :x=>X(rabif)], Ψ0=gs, mint=0.0, maxt=max_time, numt = ev_coef * num_time, saveGraph=false, saveData=false, asymptotics=false)
    jx = ev[1, :] ./ rabif.j
    jy = ev[2, :] ./ rabif.j
    jz = ev[3, :] ./ rabif.j

    p = ev[4, :]
    q = ev[5, :]

    x = sqrt(8) * rabif.λ * q
    y = -sqrt(8) * rabif.λ * rabif.δ * p
    z = 1 / Int(2 * rabif.j)

    " Normalization "
    n = sqrt.(x .* x .+ y .* y .+ z .* z)
    x0 = x ./ n
    y0 = y ./ n
    z0 = z ./ n

    # n = sqrt.(jx .* jx .+ jy .* jy .+ jz .* jz)
    # jx = jx ./ n
    # jy = jy ./ n
    # jz = jz ./ n

    colours = palette(:auto)
    num_visible_frames = Int(show_time * num_time / max_time)

    for j = 1:num_time
        pa = PlotSphere()
    
        maxj = ev_coef * (j - 1) + 1
        minj = max(maxj - ev_coef * num_visible_frames, 1)
        
        pa = PlotLine(pa, jx[minj:maxj], jy[minj:maxj], jz[minj:maxj], color=:black, lw=1, marker=:diamond)
        pa = PlotLine(pa, x0[minj:maxj], y0[minj:maxj], z0[minj:maxj], color=:red, lw=4, marker=:square)

        minj = max(j - num_visible_frames, 1)
        maxj = j

        colourIndex = 3

        for i = 1:Int(2 * rabii.j + 1)
            q = trajectory[2*i][minj:maxj]
            p = trajectory[2*i + 1][minj:maxj]

            x = sqrt(8) * λf * q
            y = -sqrt(8) * λf * rabif.δ * p
            z = 1 / Int(2 * rabii.j)

            n = sqrt.(x .* x .+ y .* y .+ z .* z)
            x = x ./ n
            y = y ./ n
            z = z ./ n

            pa = PlotLine(pa, x, y, z, color=colours[colourIndex], lw=3)
            colourIndex += 1
        end


        pa = PlotAxis(pa)
        # display(pa)

        savefig(pa, PATH * "sphere_$j.png")
    end

end

SphereAnimation(Rabi(R=50, λ=1.5, δ=0.5, j=1//2); λf=-2.0, max_time=120)
