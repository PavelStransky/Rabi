using StatsBase
using DifferentialEquations
using LinearAlgebra

include("../../Rabi.jl")


function CalculatePoints(rabii::Rabi; λf=-0.37)	
    λi = rabii.λ

    p = 0
    q = -sqrt(0.5 * (λi^2 - 1 / λi^2))

    x0 = sqrt(8) * λi * q
    y0 = -sqrt(8) * λi * rabif.δ * p
    z0 = 1 / Int(2 * rabif.j)

    " Normalization "
    n = sqrt.(x0 .* x0 .+ y0 .* y0 .+ z0 .* z0)
    x0 = x0 ./ n
    y0 = y0 ./ n
    z0 = z0 ./ n

    x = sqrt(8) * λf * q
    y = -sqrt(8) * λf * rabif.δ * p
    z = 1 / Int(2 * rabif.j)

    " Normalization "
    n = sqrt.(x .* x .+ y .* y .+ z .* z)
    x = x ./ n
    y = y ./ n
    z = z ./ n

    println(1.2 * [x0, y0, z0])
    println(1.2 * [x, y, z])
end

CalculatePoints(Rabi(R=30, λ=1.5, δ=0.5, j=1//2); λf=-sqrt(2)/5)
