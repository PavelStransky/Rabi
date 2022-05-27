@testset "Wigner" begin    
    @test begin
        rabi = Rabi(N=500, R=100, λ=0.75, δ=0.5)
        Wigner(rabi; ts=LinRange(0, 10, 11), xs=LinRange(-1.2, 1.2, 101), ys=LinRange(-0.6, 0.6, 101), saveData=false, saveGraph=false, clims=(-0.2, 0.2))
        Wigner(rabi; ts=LinRange(0, 10, 11), xs=LinRange(-1.2, 1.2, 101), ys=LinRange(-0.6, 0.6, 101), operators=[:P=>PGS(rabi)], saveData=false, saveGraph=false, clims=(-0.2, 0.2))
    end

    @test begin
        rabi = Rabi(N=500, R=100, λ=0.75, δ=0.5, μ=0.4, ν=0.4)
        Wigner(rabi; ts=LinRange(0, 10, 11), operators=[:P_log=>PGS(rabi), :q=>X(rabi), :p=>P(rabi)], xs=LinRange(-2, 1, 101), ys=LinRange(-1, 1, 101), saveData=false, saveGraph=false, clims=(-0.2, 0.2))
    end
    
    @test begin
        rabi = Rabi(N=500, R=100, λ=0.75)
        energies, vectors = eigenstates(rabi)
        a, b = ProjectParity(rabi, vectors[1], vectors[2])
        pa = projector(a)
    
        Wigner(Rabi(rabi; λ=rabi.λ / 2); Ψ0=a, operators=[:P_r=>pa], ts=LinRange(0, 10, 11), clim=(-0.1, 0.1), saveData=false, saveGraph=false)
    end
end