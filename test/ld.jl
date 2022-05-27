@testset "Level Dynamics" begin
    @test begin 
        rabi = Rabi(R=10, λ=0.75, δ=0.5)
        LevelDynamics(rabi; ylims=(-2, 2))[end,end] ≈ 5.976640557637266
    end

    @test begin 
        rabi = Rabi(R=10, λ=1.5, δ=0.5, μ=0.4, ν=0.4)
        LevelDynamics(rabi; type=:δ, ps=LinRange(-1, 1, 201), ylims=(-2, 2))[end,end] ≈ 5.8159507817005345
    end
end