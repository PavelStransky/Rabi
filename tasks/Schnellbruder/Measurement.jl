using Optim
using Measures
using ProgressMeter

include("../../Rabi.jl")
include("ClassicalTrajectory.jl")

gr()
default(size=(1000,1000))

const PATH = "c:/results/rabi/Schnellbruder/measurement/"

"""
    Converts Cartesian coordinates of a 3D vector v to spherical coordinates.
"""
function CartesianToSpherical(v)
    x, y, z = v
    r = sqrt(x^2 + y^2 + z^2)
    θ = acos(clamp(float(z / r), -1.0, 1.0))  # polar angle
    ϕ = atan(float(y), float(x))              # azimuthal angle
    return r, θ, ϕ
end


"""
Projector onto the spin-j eigenstate with magnetic number m along direction n.

- j can be Int (e.g. 1,2,3,...) or Rational (e.g. 3//2, 5//2, ...)
- m must be one of -j, -j+1, ..., j
- θ, ϕ are the polar and azimuthal angles defining the direction n on the Bloch sphere
"""
function SpinProjector(j, θ, ϕ, m)
    b = QuantumOptics.SpinBasis(j)

    Jy = sigmay(b) / 2
    Jz = sigmaz(b) / 2

    # rotation operator R = exp(-i ϕ Jz) * exp(-i θ Jy)
    R = exp(-1im * ϕ * Jz) * exp(-1im * θ * Jy)

    # find the basis index corresponding to m (robust for integer/half-integer j)
    ms = collect(-j:1:j)                 # -j, -j+1, ..., j
    idx = findfirst(==(m), ms)
    idx === nothing && error("Invalid m=$m for j=$j. Allowed m are: $ms")

    ket_mz = basisstate(b, idx)
    ket_mn = R * ket_mz

    return projector(ket_mn)
end


"""
    Bundles together the fixed parameters shared across measurement functions:
    the quantum system, spin projection quantum number, trajectory brother indices,
    the Wigner function grid, and precomputed operators.
"""
struct MeasurementConfig{R, T, X<:AbstractRange, FockOp, SpinOp, SpinKet}
    rabif::R                # final Rabi system
    m::T                    # spin projection quantum number
    brother1::Int           # index of first trajectory brother
    brother2::Int           # index of second trajectory brother
    xs::X                   # unscaled Wigner grid in x (for plotting)
    ys::X                   # unscaled Wigner grid in y (for plotting)
    xss::X                  # scaled Wigner grid in x (for computation)
    yss::X                  # scaled Wigner grid in y (for computation)
    Id::FockOp              # identity on FockBasis(rabif)
    Jy::SpinOp              # spin Jy on SpinBasis(rabif.j)
    Jz::SpinOp              # spin Jz on SpinBasis(rabif.j)
    ket_mz::SpinKet         # spin eigenstate |m⟩_z
end

function MeasurementConfig(rabif, m; wigner_lims=1.5, wigner_mesh=101,
                           brother1=1, brother2=Int(2 * rabif.j) + 1)
    xs  = LinRange(-wigner_lims, wigner_lims, wigner_mesh)
    ys  = LinRange(-wigner_lims, wigner_lims, wigner_mesh)
    xss = xs * sqrt(Size(rabif))
    yss = ys * sqrt(Size(rabif))

    Id = one(FockBasis(rabif))

    b  = QuantumOptics.SpinBasis(rabif.j)
    Jy = sigmay(b) / 2
    Jz = sigmaz(b) / 2

    ms  = collect(-rabif.j:1:rabif.j)
    idx = findfirst(==(m), ms)
    idx === nothing && error("Invalid m=$m for j=$(rabif.j). Allowed values: $ms")
    ket_mz = basisstate(b, idx)

    return MeasurementConfig(rabif, m, brother1, brother2, xs, ys, xss, yss, Id, Jy, Jz, ket_mz)
end


"""
    Projector onto the spin-m eigenstate along direction (θ, ϕ), using precomputed operators from cfg.
"""
function SpinProjector(cfg::MeasurementConfig, θ, ϕ)
    R = exp(-1im * ϕ * cfg.Jz) * exp(-1im * θ * cfg.Jy)
    return projector(R * cfg.ket_mz)
end


"""
    Calculates the direction orthogonal to the two brothers defined in cfg.
    point: classical trajectory snapshot (first element is time, then pairs of q,p per brother)
"""
function OrthogonalDirection(cfg::MeasurementConfig, point)
    rabif = cfg.rabif

    q = point[2 * cfg.brother1]
    p = point[2 * cfg.brother1 + 1]
    b1 = [rabif.λ / sqrt(2) * q, -rabif.λ * rabif.δ / sqrt(2) * p, 0.5]

    q = point[2 * cfg.brother2]
    p = point[2 * cfg.brother2 + 1]
    b2 = [rabif.λ / sqrt(2) * q, -rabif.λ * rabif.δ / sqrt(2) * p, 0.5]

    # Vector orthogonal to the two brothers (b1 x b2)
    return [b1[2] * b2[3] - b1[3] * b2[2], b1[3] * b2[1] - b1[1] * b2[3], b1[1] * b2[2] - b1[2] * b2[1]]
end


"""
    Measures the spin projection in direction given by angles θ and ϕ,
    and returns the Wigner negativity, probability of the projection, and the Wigner function.
"""
function WignerNegativity(cfg::MeasurementConfig, Ψ, θ, ϕ)
    Pr = tensor(cfg.Id, SpinProjector(cfg, θ, ϕ))

    p  = real(expect(Pr, Ψ))
    Φ  = Pr * Ψ
    Φ  = Φ / norm(Φ)
    ρq = ptrace(Φ, 2)
    w  = wigner(ρq, cfg.xss, cfg.yss)

    return 0.5 * (sum(abs.(w)) / sum(w) - 1), p, w
end


"""
    Computes the Wigner negativity before measurement, after projecting onto the orthogonal
    direction, and after optimizing over all directions.
    Returns (i, negativity, negativity_orthogonal, projection_orthogonal, negativity_maximum, projection_maximum).
"""
function Projection(cfg::MeasurementConfig, i, Ψ, point)
    if i == 1
        return i, 0.0, 0.0, 0.0, 0.0, 0.0
    end

    # Wigner negativity before measurement
    ρq = ptrace(Ψ, 2)
    w  = wigner(ρq, cfg.xss, cfg.yss)
    negativity = 0.5 * (sum(abs.(w)) / sum(w) - 1)

    # Orthogonal direction
    bc = OrthogonalDirection(cfg, point)
    bc ./= norm(bc)
    r, θ, ϕ = CartesianToSpherical(bc)

    negativity_orthogonal, projection_orthogonal, _ = WignerNegativity(cfg, Ψ, θ, ϕ)

    # Maximize negativity — start from the orthogonal direction (physically motivated, faster convergence)
    result = optimize(x -> -WignerNegativity(cfg, Ψ, x[1], x[2])[1], [θ, ϕ], NelderMead())
    negativity_maximum, projection_maximum, _ = WignerNegativity(cfg, Ψ, result.minimizer[1], result.minimizer[2])

    return i, negativity, negativity_orthogonal, projection_orthogonal, negativity_maximum, projection_maximum
end


"""
Measures the projection onto the state defined by the classical trajectory,
and the Wigner negativity before and after the measurement.
"""
function ProjectionTime(rabii, λf, m; maxt=200, numt=2001, wigner_mesh=101, wigner_lims=1.5)
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

    @showprogress for (i, Ψ) in enumerate(Ψt)
        point = getindex.(trajectory, i)
        _, negativity[i], negativity_orthogonal[i], projection_orthogonal[i], negativity_maximum[i], projection_maximum[i] = Projection(cfg, i, Ψ, point)
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


function ScanNegativity(rabii, λf; t=31.4, wigner_mesh=101, wigner_mesh_final=201, wigner_lims=1.5, scan_mesh=101, m=-rabii.j, brother1=1, brother2=Int(2 * rabii.j) + 1)
    rabif = Copy(rabii; λ=λf)
    println(rabii.N)

    Ψ0 = SingleWellState(rabii)
    ts  = LinRange(0.0, t, 2)

    negativity = Matrix{Float64}(undef, scan_mesh, scan_mesh)
    projection = Matrix{Float64}(undef, scan_mesh, scan_mesh)

    print("Starting time evolution...")
    time = @elapsed tout, Ψt = timeevolution.schroedinger(ts, Ψ0, H(rabif))
    println(time, "s")

    Ψ = Ψt[end]
    println("Final Schrödinger time (check): ", tout[end])

    cfg = MeasurementConfig(rabif, m; wigner_lims=wigner_lims, wigner_mesh=wigner_mesh, brother1=brother1, brother2=brother2)

    # Scan over all measurement directions
    ϕs = LinRange(0, 2π, scan_mesh)
    θs = LinRange(0, π, scan_mesh)

    @showprogress for i in eachindex(ϕs)
        for j in eachindex(θs)
            n, p, _ = WignerNegativity(cfg, Ψ, θs[j], ϕs[i])
            negativity[i, j] = n
            projection[i, j] = p
        end
    end

    maxidx = argmax(negativity)
    maxn = negativity[maxidx]
    maxθ = θs[maxidx[2]]
    maxϕ = ϕs[maxidx[1]]

    # Orthogonal direction from classical trajectory
    trajectory = last.(Trajectory(rabif, rabii.λ, t, 2))
    println("Trajectory time (check): ", trajectory[1])
    bc = OrthogonalDirection(cfg, trajectory)
    bc ./= norm(bc)
    r, θo, ϕo = CartesianToSpherical(bc)

    # Plot scan results
    p = contourf(ϕs, θs, transpose(negativity), title="Negativity $(t) $brother1$brother2", c=:turbo, grid=false, xlabel="ϕ", ylabel="θ", left_margin=10mm, bottom_margin=5mm)
    p = scatter!(p, [ϕo], [θo], markershape=:star5, markersize=15, label="Orthogonal direction", color=:white)
    p = scatter!(p, [maxϕ], [maxθ], markershape=:circle, markersize=15, label="Max negativity $(maxn)", color=:white)
    display(p)
    savefig(p, "$(PATH)negativity_scan_$(rabii)_$(λf)_$(t)_$(Int(m + rabii.j))_$brother1$brother2.png")

    p = contourf(ϕs, θs, transpose(projection), title="Projection $(t) $brother1$brother2", c=:turbo, grid=false, xlabel="ϕ", ylabel="θ", left_margin=10mm, bottom_margin=5mm)
    p = scatter!(p, [ϕo], [θo], markershape=:star5, markersize=15, label="Orthogonal direction", color=:white)
    p = scatter!(p, [maxϕ], [maxθ], markershape=:circle, markersize=15, label="Max negativity $(maxn)", color=:white)
    display(p)
    savefig(p, "$(PATH)projection_scan_$(rabii)_$(λf)_$(t)_$(Int(m + rabii.j))_$brother1$brother2.png")

    # Final Wigner plots at higher resolution
    cfg_final = MeasurementConfig(rabif, m; wigner_lims=wigner_lims, wigner_mesh=wigner_mesh_final, brother1=brother1, brother2=brother2)

    maxn, maxp, w = WignerNegativity(cfg_final, Ψ, maxθ, maxϕ)
    clim = (-maximum(abs.(w)), maximum(abs.(w)))
    p = heatmap(cfg_final.xs, cfg_final.ys, transpose(w), title="Maximum $(t) $brother1$brother2", c=:bwr, grid=false, xlabel="x", ylabel="p", left_margin=10mm, bottom_margin=5mm, clim=clim)
    display(p)
    savefig(p, "$(PATH)wigner_max_$(rabii)_$(λf)_$(t)_$(Int(m + rabii.j))_$brother1$brother2.png")

    no, po, w = WignerNegativity(cfg_final, Ψ, θo, ϕo)
    clim = (-maximum(abs.(w)), maximum(abs.(w)))
    p = heatmap(cfg_final.xs, cfg_final.ys, transpose(w), title="Orthogonal $(t) $brother1$brother2", c=:bwr, grid=false, xlabel="x", ylabel="p", left_margin=10mm, bottom_margin=5mm, clim=clim)
    display(p)
    savefig(p, "$(PATH)wigner_orthogonal_$(rabii)_$(λf)_$(t)_$(Int(m + rabii.j))_$brother1$brother2.png")

    Export("$(PATH)negativity_scan_$(rabii)_$(λf)_$(t)_$(Int(m + rabii.j))_$brother1$brother2", negativity)
    Export("$(PATH)projection_scan_$(rabii)_$(λf)_$(t)_$(Int(m + rabii.j))_$brother1$brother2", projection)
    Export("$(PATH)wigner_orthogonal_$(rabii)_$(λf)_$(t)_$(Int(m + rabii.j))_$brother1$brother2", cfg_final.xs, cfg_final.ys, w)
    Export("$(PATH)wigner_max_$(rabii)_$(λf)_$(t)_$(Int(m + rabii.j))_$brother1$brother2", cfg_final.xs, cfg_final.ys, w)

    return maxθ, maxϕ, maxn, maxp, ϕo, θo, no, po
end


function ScanNegativityTime(rabii, λf; maxt=100.0, numt=101, wigner_mesh=101, scan_mesh=101)
    ts = LinRange(0.0, maxt, numt)

    maxθs = Vector{Float64}(undef, numt)
    maxϕs = Vector{Float64}(undef, numt)
    maxns = Vector{Float64}(undef, numt)
    maxps = Vector{Float64}(undef, numt)
    ϕos   = Vector{Float64}(undef, numt)
    θos   = Vector{Float64}(undef, numt)
    nos   = Vector{Float64}(undef, numt)
    pos   = Vector{Float64}(undef, numt)

    for (i, t) in enumerate(ts)
        if t == 0
            maxθs[i] = maxϕs[i] = maxns[i] = maxps[i] = 0.0
            ϕos[i]   = θos[i]   = nos[i]   = pos[i]   = 0.0
            continue
        end
        println("Time: ", t)
        maxθs[i], maxϕs[i], maxns[i], maxps[i], ϕos[i], θos[i], nos[i], pos[i] = ScanNegativity(rabii, λf; t=t, wigner_mesh=wigner_mesh, scan_mesh=scan_mesh, m=-rabii.j)
    end

    p = plot(ts, maxns, xlabel="\$t\$", ylabel="Max Negativity", title="Max Wigner Negativity over time", legend=false)
    display(p)

    p = plot(ts, maxθs, xlabel="\$t\$", label="Max θ", legend=true)
    p = plot!(p, ts, maxϕs, xlabel="\$t\$", label="Max ϕ")
    p = plot!(p, ts, θos, xlabel="\$t\$", label="Orthogonal θ")
    p = plot!(p, ts, ϕos, xlabel="\$t\$", label="Orthogonal ϕ")
    display(p)

    p = plot(ts, maxps, xlabel="\$t\$", label="Max norm")
    p = plot!(p, ts, pos, xlabel="\$t\$", label="Orthogonal norm")
    display(p)

    Export("$(PATH)maxtheta_$(rabii)_$(λf)", maxθs)
    Export("$(PATH)maxphi_$(rabii)_$(λf)", maxϕs)
    Export("$(PATH)maxns_$(rabii)_$(λf)", maxns)
    Export("$(PATH)maxps_$(rabii)_$(λf)", maxps)
    Export("$(PATH)no_$(rabii)_$(λf)", nos)
    Export("$(PATH)phio_$(rabii)_$(λf)", ϕos)
    Export("$(PATH)thetao_$(rabii)_$(λf)", θos)
    Export("$(PATH)pos_$(rabii)_$(λf)", pos)
end

# t = Projection(Rabi(R=20, λ=1.5, δ=0.5, j=1//2), -0.283, -1//2, maxt=100, numt=1001, wigner_mesh=201)
# t = Measure(Rabi(R=20, λ=1.5, δ=0.5, j=4//2), -0.283, -4//2, maxt=100, numt=1001, wigner_mesh=201)
# MaximizeNegativityDirection(Rabi(R=20, λ=1.5, δ=0.5, j=1//2), -0.283, t=31.4, num=200, mesh=401)
# MaximizeNegativityDirection(Rabi(R=20, λ=1.5, δ=0.5, j=4//2), -0.283, t=47.1, num=1000, mesh=401)
# ScanNegativity(Rabi(R=20, λ=1.5, δ=0.5, j=4//2), -0.283, m=4//2, t=31.4, wigner_mesh=101, scan_mesh=51)

# ScanNegativityTime(Rabi(R=20, λ=1.5, δ=0.5, j=1//2), -0.283, numt=201, wigner_mesh=101, scan_mesh=101)
# ScanNegativityTime(Rabi(R=20, λ=1.5, δ=0.5, j=1//2), -0.283, numt=6, wigner_mesh=31, scan_mesh=31)

# for m in -4//2:1:4//2
#     ScanNegativity(Rabi(R=20, λ=1.5, δ=0.5, j=4//2), -0.283, m=m, t=47.1, wigner_mesh=101, scan_mesh=101, brother1=1, brother2=5)
# end

# t = Projection(Rabi(R=20, λ=1.5, δ=0.5, j=4//2), -0.283, -4//2, maxt=200, numt=1001, wigner_mesh=151)
# t = Projection(Rabi(R=20, λ=1.5, δ=0.5, j=4//2), -0.283, -2//2, maxt=200, numt=1001, wigner_mesh=151)
# t = Projection(Rabi(R=20, λ=1.5, δ=0.5, j=4//2), -0.283, -0//2, maxt=200, numt=1001, wigner_mesh=151)

# print(ScanNegativity(Rabi(R=20, λ=1.5, δ=0.5, j=4//2), -0.283, m=-4//2, t=31.4, wigner_mesh=101, scan_mesh=11, brother1=2, brother2=4, save=false))

# Supplemental material
# Projection(Rabi(R=50, λ=1.5, δ=0.5, j=1//2), -0.369, -1//2, maxt=200, numt=1001, wigner_mesh=201)
# Projection(Rabi(R=50, λ=1.5, δ=0.5, j=2//2), -0.369, -2//2, maxt=200, numt=1001, wigner_mesh=201)

# Projection(Rabi(R=20, λ=1.5, δ=0.5, j=2//2), -0.369, -2//2, maxt=200, numt=101, wigner_mesh=51)
