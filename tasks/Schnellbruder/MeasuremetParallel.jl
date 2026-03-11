using Distributed

include("../../Rabi.jl")
include("ClassicalTrajectory.jl")

WORKERS = 6

if nprocs() <= WORKERS
    addprocs(WORKERS + 1 - nprocs())
end

@everywhere include("Measurement.jl")

"""
Measures the projection onto the state defined by the classical trajectory,
and the Wigner negativity before and after the measurement — distributed across workers.
"""
function ProjectionParallel(rabii, λf, m; maxt=200, numt=2001, wigner_mesh=101, wigner_lims=1.5)
    rabif = Copy(rabii; λ=λf)
    println(rabii.N)

    Ψ0 = SingleWellState(rabii)
    ts  = LinRange(0.0, maxt, numt)

    negativity            = Vector{Float64}(undef, numt)
    negativity_orthogonal = Vector{Float64}(undef, numt)
    projection_orthogonal = Vector{Float64}(undef, numt)
    negativity_maximum    = Vector{Float64}(undef, numt)
    projection_maximum    = Vector{Float64}(undef, numt)

    print("Starting time evolution...")
    time = @elapsed tout, Ψt = timeevolution.schroedinger(ts, Ψ0, H(rabif))
    println(time, "s")

    cfg        = MeasurementConfig(rabif, m; wigner_lims=wigner_lims, wigner_mesh=wigner_mesh)
    trajectory = Trajectory(rabif, rabii.λ, maxt, numt)

    input  = [(i, Ψt[i], getindex.(trajectory, i)) for i in eachindex(Ψt)]
    result = progress_pmap(x -> Projection(cfg, x[1], x[2], x[3]), input)

    for i in eachindex(result)
        j = result[i][1]

        negativity[j]            = result[i][2]
        negativity_orthogonal[j] = result[i][3]
        projection_orthogonal[j] = result[i][4]
        negativity_maximum[j]    = result[i][5]
        projection_maximum[j]    = result[i][6]
    end

    pl = plot(ts, negativity, xlabel="\$t\$", ylabel="Original", title="Evolution of Wigner negativity", legend=true)
    pl = plot!(pl, ts, negativity_orthogonal, xlabel="\$t\$", ylabel="Measured")
    pl = plot!(pl, ts, negativity_maximum, xlabel="\$t\$", ylabel="Maximum")
    display(pl)

    p2 = plot(ts, projection_orthogonal, xlabel="\$t\$", ylabel="P", title="Orthogonal")
    p2 = plot!(p2, ts, projection_maximum, xlabel="\$t\$", label="Maximum")
    display(p2)

    Export("$(PATH)negativity_$(rabii)_$(λf)_$(Int(m + rabii.j))", tout, negativity)
    Export("$(PATH)negativity_orthogonal_$(rabii)_$(λf)_$(Int(m + rabii.j))", tout, negativity_orthogonal)
    Export("$(PATH)negativity_maximum_$(rabii)_$(λf)_$(Int(m + rabii.j))", tout, negativity_maximum)
    Export("$(PATH)projection_orthogonal_$(rabii)_$(λf)_$(Int(m + rabii.j))", tout, projection_orthogonal)
    Export("$(PATH)projection_maximum_$(rabii)_$(λf)_$(Int(m + rabii.j))", tout, projection_maximum)
end

# ProjectionParallel(Rabi(R=20, λ=1.5, δ=0.5, j=1//2), -0.369, -1//2, maxt=100, numt=2001, wigner_mesh=101)
ProjectionParallel(Rabi(R=20, λ=1.5, δ=0.5, j=2//2), -0.369, -2//2, maxt=100, numt=2001, wigner_mesh=101)