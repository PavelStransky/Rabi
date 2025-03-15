include("../Rabi.jl")

const PATH = "d:/results/rabi/wf/"

pyplot(size=(600, 600))


""" Evolution of the partial trace """
function PartialTraceEvolution(mint=0.0, maxt=100.0, numt=1000)
    R = 100
    # Initial state

    λf = 2.5

    # Initial and final systems
    rabif = Rabi(R=R, λ=λf)
    rabif = Rabi(R=R, λ=λf, δ=0.5, μ=0.001, ν=0.001)
    gs = ΨGS(rabif)
   
    spin = Array{Float64}(undef, numt, 4)
    ptt = Array{Float64}(undef, numt, 4)
    pt = Array{Float64}(undef, numt, 2)
    field = Array{Float64}(undef, numt, 4)

    time = @elapsed begin
        u = U(rabif, mint)
        du = U(rabif, (maxt - mint) / numt)
        Ψ = u * gs      # Initial state

        for i in 1:numt
            ρq = ptrace(Ψ, 1)
            println(ρq.data)

            pt[i, 1] = real(ρq.data[1, 1])
            pt[i, 2] = real(ρq.data[2, 2])

            ptt[i, 1] = real(ρq.data[1, 2])
            ptt[i, 2] = -imag(ρq.data[1, 2])
            ptt[i, 3] = real(ρq.data[1, 1]) - 0.5
            ptt[i, 4] = sqrt(ptt[i, 1]^2 + ptt[i, 2]^2 + ptt[i, 3]^2)

            q = ExpectationValue("X", X(rabif), Ψ, rabif)
            p = ExpectationValue("P", P(rabif), Ψ, rabif)

            field[i, 1] = sqrt(8) * rabif.λ * q
            field[i, 2] = -sqrt(8) * rabif.λ * rabif.δ * p
            field[i, 3] = 1 / (2 * rabif.λ) + sqrt(8) * rabif.μ * q
            field[i, 4] = sqrt(field[i, 1]^2 + field[i, 2]^2 + field[i, 3]^2)

            spin[i, 1] = ExpectationValue("Jx", Jx(rabif), Ψ, rabif)
            spin[i, 2] = ExpectationValue("Jy", Jy(rabif), Ψ, rabif)
            spin[i, 3] = ExpectationValue("Jz", Jz(rabif), Ψ, rabif)
            spin[i, 4] = sqrt(spin[i, 1]^2 + spin[i, 2]^2 + spin[i, 3]^2)

            Ψ = du * Ψ
        end
    end
    println(time)

    tout = LinRange(mint, maxt, numt)

    p = plot(tout, ptt[:, 1:3], xlabel="\$t\$", title="Partial trace", legend=false)
    display(p)
    savefig(p, "$(PATH)ptt_($(rabif)).png")

    p = plot(tout, pt, xlabel="\$t\$", title="Partial trace - diagonal elements", legend=false)
    display(p)
    savefig(p, "$(PATH)pt_($(rabif)).png")

    p = plot(tout, field[:, 1:3], xlabel="\$t\$", title="B", legend=false)
    display(p)
    savefig(p, "$(PATH)field_($(rabif)).png")

    p = plot(tout, spin[:, 1:3], xlabel="\$t\$", title="J", legend=false)
    display(p)
    savefig(p, "$(PATH)spin_($(rabif)).png")

    println(ptt[:, 4])
end

PartialTraceEvolution()

# LD(4//2)

# for j in 1:20
#     StrengthFunction(j // 2)
# end

# PartialTrace(4//2)