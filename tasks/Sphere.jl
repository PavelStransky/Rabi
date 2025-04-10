using StatsBase
using DifferentialEquations
using LinearAlgebra

# const PATH = "d:/results/rabi/schnellbruder/"
const PATH = "/home/stransky/results/"

gr()
default(size=(600,600), dpi=300)

function InitialState(rabii)
    _, vs = eigenstates(rabii, 2)    
    a, b = ProjectParity(rabii, vs[1], vs[2])
    gs1 = (a + b) / sqrt(2)
    gs2 = (a - b) / sqrt(2)

    q1 = ExpectationValue("E", X(rabii), gs1, rabii)
    q2 = ExpectationValue("E", X(rabii), gs2, rabii)

    if q1 < 0
        return gs1
    end

    return gs2
end

function EquationOfMotion!(dx, x, parameters, t)
    q, p = x
    rabi, m = parameters

    s2 = 2 * rabi.λ^2 * (q^2 + rabi.δ^2 * p^2) + 1
    if s2 < 0                       
        s2 = 0                      # With isoutofdomain=CheckDomain the calculation shouldn't enter here, but in enters anyway
        @error("s negative!") 
    end
    
    s = sqrt(s2)

    dx[1] = p * (1 + rabi.λ^2 * rabi.δ^2 * m / rabi.j / s)
    dx[2] = -q * (1 + rabi.λ^2 * m / rabi.j / s)
end

function Trajectory(rabif, t, λi)
    x0 = [-sqrt(0.5 * (λi^2 - 1 / λi^2)), 0]
    timeInterval = (0.0, t)
    solver = TsitPap8()
    tolerance = 1E-6
    fnc = ODEFunction(EquationOfMotion!)

    result = []

    for m = -rabif.j:rabif.j
        problem = ODEProblem(fnc, x0, timeInterval, (rabif, m))

        saveat = collect(range(max(t - 10, 0), t, step=0.05))
        time = @elapsed solution = solve(problem, solver, reltol=tolerance, abstol=tolerance, verbose=true, saveat=saveat)

        push!(result, collect[solution])
    end

    return p
end


function VisibilityOnSphere(point)
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

    p = surface(x, y, z, alpha=0.1, legend=false, color=:viridis, camera=(30,30), 
    grid=false, ticks=nothing, axis=false,
    xlims=(-lim, lim), ylims=(-lim, lim), zlims=(-lim, lim))

    return p
end

function PlotLine(p, xs, ys, zs; color=:red, lw=2)
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
            p = plot!(p, lx, ly, lz, color=:red, alpha=if visible 1.0 else 0.3 end, lw=lw)
            
            visible = v
            lx = [x]
            ly = [y]
            lz = [z]
        end
    end

    p = plot!(p, lx, ly, lz, color=:red, alpha=if visible 1.0 else 0.3 end, lw=lw)

    return p
end

function pokus(rabif, t, p, λi)
    x = []
    y = []
    z = []

    for ph = range(0, stop=2π, length=50)
        # Define a point on the sphere (example: north pole)
        th = 0.7
        px, py, pz = sin(th) * cos(ph), sin(th) * sin(ph), cos(th)

        push!(x, px)
        push!(y, py)
        push!(z, pz)
    end

    pa = PlotSphere()
    pa = PlotAxis(pa)
    pa = PlotLine(pa, x, y, z, color=:red, lw=4)
    # pa = plot(pa, margin=-5mm)

    p = plot(p, pa, layout=grid(2, 1, heights=[0.3,0.7]))

    return p
end

function pokus()
    λf = -0.37
    rabii = Rabi(R=50, λ=1.5, δ=0.5, j=2//2)
    rabif = Copy(rabii, λ=λf)

    return Trajectory(rabif, 10, rabii.λ)
end

print(pokus())



