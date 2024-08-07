include("Export.jl")
include("Rabi.jl")
include("QuarticOscillator.jl")

pyplot(size=(1000, 1000))

const PATH = "d:\\results\\Rabi\\oscillator\\"

function RabiCat()
    R = 100
    λ = 0.75
    μ = 0.001
    ν = μ
    n = 1000
    num = 401

    rabi = Rabi(N=n, R=R, λ=λ, μ=μ, ν=ν)

    ts = LinRange(0, 100, 2001)

    #Ψ0 = tensor(squeeze(FockBasis(rabi), 2) * coherentstate(FockBasis(rabi), 0), spindown(SpinBasis(rabi)))

    rabi0 = Rabi(N=n, R=R)
    Ψ0 = eigenstates(rabi0)
    Ψ0 = Ψ0[2][2]

    result = ExpectationValues(rabi, AllOperators(rabi); Ψ0=Ψ0, mint=0.0, maxt=50.0, numt=1000, showGraph=true, saveGraph=true, saveData=true, asymptotics=false)

    sf = StrengthFunction(rabi, Ψ0)

    p = scatter(sf[1], sf[2], markeralpha=0.7, markerstrokewidth=0, markersize=6, xlims=(-2, 2))
    display(p)

    #_, p = Overlap(vectors; limit=length(vectors) ÷ 3)
    #savefig(p, "$(PATH)psi_($(rabi.R),$(rabi.μ)).png")

    # Level Dynamics
    #_, pld = LevelDynamics(rabi, ps=LinRange(0, 1.2, 601), limit=300, ylims=(-1,1))

    #p = scatter(pld, λ .- 1.0 .* sf[2], sf[1], markeralpha=0.9, markerstrokewidth=0, markersize=8)
    #p = plot(p, [0, λ], [-0.5, ExpectationValue(H(rabi), Ψ0) / rabi.R], linewidth=3, arrow=(:closed, 1.0))
    #display(plot(p))
end

function CUSPCat()
    for a = LinRange(0, 0.02, 101)
        cusp = CUSP(N=2000, ħ=0.1, a=a)
    #_, pld = LevelDynamics(cusp, ps=LinRange(-4, 1, 401), type=:b, limit=100, ylims=(-2,2))

        ts = LinRange(0, 1500, 301)

        Ψ0 = coherentstate(Basis(cusp), 0)
        result = ExpectationValues(cusp, AllOperators(cusp); Ψ0=Ψ0, mint=0.0, maxt=1500.0, numt=301, showGraph=true, saveGraph=true, saveData=true, asymptotics=false)
    end
    
    #sf = StrengthFunction(cusp, Ψ0)
    #p = scatter(sf[1], sf[2], markeralpha=0.7, markerstrokewidth=0, markersize=6, xlims=(-2, 2))
    #display(p)
end

CUSPCat()
#RabiCat()