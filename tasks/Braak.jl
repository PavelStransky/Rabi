using Plots

""" Calculations based on Braak2011 paper """

include("../Rabi.jl")

PATH = "d:/results/rabi/braak/"

g(rabi::Rabi) = rabi.λ * sqrt(rabi.R)
Δ(rabi::Rabi) = 0.5 * rabi.ω * rabi.R

f(rabi::Rabi, n, x) = 2.0 * g(rabi) / rabi.ω + 0.5 / g(rabi) * (n * rabi.ω - x + Δ(rabi)^2 / (x - n * rabi.ω))
fcache(rabi::Rabi, maxn, x) = [f(rabi, n, x) for n in 0:maxn]

function K(rabi::Rabi, n, x)
    if n == 0
        return 1
    elseif n == 1
        return f(rabi, 0, x)
    else
        return (f(rabi, n - 1, x) * K(rabi, n - 1, x) - K(rabi, n - 2, x)) / n
    end
end
function Kcache(fcache)
    nmax = length(fcache)
    result = Array{Float64}(undef, nmax + 1)
    result[1] = 1
    result[2] = fcache[1]
    
    for n in 2:nmax
        result[n + 1] = (fcache[n] * result[n] - result[n - 1]) / n
    end

    return result
end

function G(rabi::Rabi, sign, x; cutoff=30)
    fc = fcache(rabi, cutoff, x)
    Kc = Kcache(fc)

    result = 0.0
    for n in 0:cutoff
        result += Kc[n + 1] * (1 - sign * Δ(rabi) / (x - n * rabi.ω)) * (g(rabi) / rabi.ω)^n
    end

    if abs(result) > 1E3 result = NaN end
    return result
end


xs = LinRange(-10, 20, 120000)
limit = 1
R = 2.8

λs = LinRange(0, 1.5, 501)
ld, pld = LevelDynamics(Rabi(R=R), ps=λs, limit=100, ylims=(-3,2))

pyplot(size=(1920,1080))

for (j, λ) in enumerate(λs)
    rabi = Rabi(R=R, λ=λ)

    xss = xs .- (g(rabi)^2 / rabi.ω)
    xss = xss ./ (rabi.R)

    p1 = plot([-10, 5], [0, 0], color=:gray)
    p1 = plot!(p1, xss, [G(rabi, 1.0, x) for x in xs], lw=2, ylims=(-limit, limit), xlims=(-3, 2))
    p1 = plot!(p1, xss, [G(rabi, -1.0, x) for x in xs], lw=2, legend=false)

    p2 = plot(pld)
    p2 = plot!(p2, [λ, λ], [-3, 2])
    p2 = scatter!(p2, LinRange(λ, λ, length(ld[j, :])), ld[j,:])

    p=plot(p2, p1)
    display(p)
    savefig(p, "$(PATH)braak$j.png")
end