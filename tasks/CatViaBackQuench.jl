include("../Rabi.jl")

const PATH = "d:/results/rabi/wf/"

pyplot(size=(600, 600))

function LD(j)
    R = 10
    _, pld = LevelDynamics(Rabi(R=R, N=300, j=j), ps=LinRange(0, 2.5, 501), limit=500, ylims=(-1,2), saveGraph=true)
end

function StrengthFunction(j)
    R = 50

    λi = 1.5
    λf = 0.75

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

for j in 1:20
    CalculateJ(j // 2)
end
