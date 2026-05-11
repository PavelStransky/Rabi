include("../../Rabi.jl")

function ProjectParity(rabi::Rabi, vector1::Ket, vector2::Ket)
    p = Parity(rabi)
    m = real([dagger(vector1) * p * vector1 dagger(vector1) * p * vector2;
        dagger(vector2) * p * vector1 dagger(vector2) * p * vector2])
    
    ev = eigen(m)

    a = vector1 * ev.vectors[1, 2] + vector2 * ev.vectors[2, 2]
    b = vector1 * ev.vectors[1, 1] + vector2 * ev.vectors[2, 1]

    return a, b

end

" A parity-violating state for a backward quench "
function SingleWellState(rabi::Rabi; left=true)
    _, vs = eigenstates(rabi, 2) 
    a, b = ProjectParity(rabi, vs[1], vs[2])

    gs1 = (a + b) / sqrt(2)
    gs2 = (a - b) / sqrt(2)

    q1 = ExpectationValue("E", X(rabi), gs1, rabi)

    if q1 < 0 && left
        return gs1
    end

    return gs2
end

function CalculateWigner(R, δ, λfs; wignerMesh=401, limits=3.0)
    xs = LinRange(-limits, limits, wignerMesh)
    ys = LinRange(-1, 1, wignerMesh)

    wn = zeros(length(λfs))

    w = 0
    for (i, λf) in enumerate(λfs)
         # Initial and final systems
        rabi = Rabi(R=R, δ=δ, λ=λf, j=1//2)

        es, vs = eigenstates(rabi, 2)
        a, b = ProjectParity(rabi, vs[1], vs[2])
      
        xss = xs * sqrt(Size(rabi))
        yss = ys * sqrt(Size(rabi))

        ρ = DensityMatrix(rabi, a)
        w = wigner(ρ, xss, yss)

        wn[i] = 0.5 * (sum(abs, w) * (xss[2] - xss[1]) * (yss[2] - yss[1]) - 1)
        println("λf: $λf, Wigner Negativity: $(wn[i])")
    end

    p = plot(xs, ys, w', st=:heatmap, title="Wigner Function", xlabel="x", ylabel="y")
    display(p)

    p = plot(λfs, wn, label="Wigner Negativity", xlabel="λf", ylabel="Wigner Negativity", title="Wigner Negativity vs λf")
    display(p)

    return wn
end

wn = CalculateWigner(20, 0.0, LinRange(0, 3, 121))
