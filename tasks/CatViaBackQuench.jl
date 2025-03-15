include("../Rabi.jl")

const PATH = "d:/results/rabi/wf/"

pyplot(size=(600, 600))

function LD(j)
    R = 10
    _, pld = LevelDynamics(Rabi(R=R, N=500, j=j), ps=LinRange(0, 2.5, 501), limit=500, ylims=(-1,2), saveGraph=true)
end

function StrengthFunction(j)
    R = 50

    λi = 1.5
    λf = 0.0

    showGraph = false
    wignerMesh = 401

    # Initial and final systems
    rabii = Rabi(R=R, λ=λi, j=j)
    rabif = Copy(rabii; λ=λf)

    ei, vi = eigenstates(rabii)
    gs = vi[1]

    # Exact parity states
    a, b = ProjectParity(rabii, vi[1], vi[2])
    c = sqrt(0.5) * (a + b)
    d = sqrt(0.5) * (a - b)

    println("GS ", ExpectationValue("parity GS", Parity(rabii), gs, rabii))
    println("A ", ExpectationValue("parity A", Parity(rabii), a, rabii))
    println("C ", ExpectationValue("parity C", Parity(rabii), c, rabii))

    # Strength function
    sfa = StrengthFunction(rabif, a)

    Export("$(PATH)sf($(Int(2*rabif.j)),$(rabif.R),$(rabif.λ),$(rabif.δ),$(rabif.μ))", sfa[1], sfa[2])
    p = scatter(sfa[1], log10.(sfa[2]),  xlims=(0,3), ylims=(-20, 0), title="$(j)", markeralpha=0.5, markerstrokewidth=0, xlabel="\$E\$", ylabel=raw"$\log_{10}S$", label=raw"$|\psi_{\mathrm{GS}}\rangle$")
    display(plot(p))

    savefig(p, "$(PATH)sf($(Int(2*rabif.j)),$(rabif.R),$(rabif.λ),$(rabif.δ),$(rabif.μ)).png")
end

function PartialTrace(j)
    R = 20

    λs = LinRange(0, 5.0, 501)
    trace = Array{Float64}(undef, 0)    
    eigen = Array{Float64}(undef, 0)

    time = @elapsed begin
        for λ in λs
            print(".")
            system = Rabi(R=R, λ=λ, j=j)
            es, vs = eigenstates(system, 2)

            a, b = ProjectParity(system, vs[1], vs[2])
            
            ρq = ptrace(a, 1)
            ev, vs = QuantumOptics.eigenstates(ρq)

            append!(trace, real(diag(ρq.data)))
            append!(eigen, ev)
        end
    end
    println(time, "s")

    trace = transpose(reshape(trace, Int(2*j+1), :))
    eigen = transpose(reshape(eigen, Int(2*j+1), :))
    p = plot(λs, trace, legend=false, color=RGBA(0, 0, 0, 0.5), xlabel=raw"$\lambda$", ylabel=raw"$\rho$")
    p = plot!(p, λs, eigen, color=RGBA(1, 0, 0, 0.5))
    display(p)

    savefig(p, "$(PATH)partialtrace_(N=$(Int(2*j)),R=$(R).png")
    println(trace[end, :])
end

""" Evolution of the partial trace """
function PartialTraceEvolution(mint=0.0, maxt=10.0, numt=1000)
    R = 50
    # Initial state

    λi = 1.5
    λf = 0.5

    # Initial and final systems
    rabii = Rabi(R=R, λ=λi, μ=0.1, ν=0.1)
    rabii = Rabi(R=R, λ=λi)
    rabif = Copy(rabii; λ=λf)

    ei, vi = eigenstates(rabii)
    gs = vi[1]

    # Exact parity states
    a, b = ProjectParity(rabii, vi[1], vi[2])
    c = sqrt(0.5) * (a + b)
    d = sqrt(0.5) * (a - b)
    
    result = Array{Float64}(undef, numt, 4)

    time = @elapsed begin
        u = U(rabif, mint)
        du = U(rabif, (maxt - mint) / numt)
        Ψ = u * a      # Initial state

        for i in 1:numt
            ρq = ptrace(Ψ, 1)
            println(ρq.data)
            result[i, 1] = 2*real(ρq.data[1, 2])
            result[i, 2] = 2*imag(ρq.data[1, 2])
            result[i, 3] = 2*real(ρq.data[1, 1]) - 1
            result[i, 4] = sqrt(result[i, 1]^2 + result[i, 2]^2 + result[i, 3]^2)
            Ψ = du * Ψ
        end
    end
    println(time)

    tout = LinRange(mint, maxt, numt)

    p = plot(tout, result, xlabel="\$t\$", title="Components of the vector", legend=false)
    display(p)

    println(result[:, 4])
end

PartialTraceEvolution()

# LD(4//2)

# for j in 1:20
#     StrengthFunction(j // 2)
# end

# PartialTrace(4//2)