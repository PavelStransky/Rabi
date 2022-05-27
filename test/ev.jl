
@testset "Expectation values" begin
    @test begin
        rabi = Rabi(R=100, λ=0.75, δ=0.5)
        ExpectationValues(rabi, AllOperators(rabi))[:,end] ≈ [0.0, 0.0, -0.28969443962354613, 0.10311964522525396, 0.0, 0.0, 61.43276562386805, 0.022718653607638927]
    end

    @test begin
        rabi = Rabi(R=100, λ=0.75, δ=0.5, μ=0.001, ν=0.001)
        ExpectationValues(rabi, AllOperators(rabi))[:,end] ≈ [0.11652134370014461, 0.049773182129196326, -0.3463066549432686, 0.030409751941296123, -0.24775081738102622, 0.18094956127002715, 36.74761215657957, 0.034929919344686484]
    end

    @test begin
        rabi = Rabi(R=100, λ=0.75, δ=0.5, μ=0.4, ν=0.4)
        ExpectationValues(rabi, AllOperators(rabi))[:,end] ≈ [0.11904081115525772, -0.005955632636221913, -0.2527271635558117, 0.1778215717195012, -0.49461480259539664, -0.16336932821414238, 93.7543218938091, 0.017269746375052325]
    end

    @test begin
        rabi = Rabi(N=500, R=100, λ=0.75)
        energies, vectors = eigenstates(rabi)
        a, b = ProjectParity(rabi, vectors[1], vectors[2])
        pa = projector(a)

        rabi = Rabi(rabi; λ=rabi.λ / 2)
        operators = AllOperators(rabi)
        operators[4] = :P => pa
        operators[8] = :r_r => pa

        ExpectationValues(rabi, operators; Ψ0=a)[:,end] ≈ [-0.08449520713403727, -0.02046440279182601, -0.4250022812926567, 0.16869206575261686, 0.2722479221300603, -0.048802142828842716, 36.89286534155856, 0.017796803222338276]
    end
end