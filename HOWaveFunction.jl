""" Logarithm of the Hermite polynomial in the form result * exp(exponent) """
function HermiteLog(n, x)
    hm1 = 2.0 * x
    hm2 = 1.0

    result = 0.0
    exponent = 0.0

    if n == 0
        result = hm2
    elseif n == 1
        result = hm1
    else
        for i = 2:n
            result = 2.0 * x * hm1 - 2 * (i - 1) * hm2
            absolute = abs(result)

            if absolute == 0
                absolute = 1
            end

            exponent += log(absolute)

            result /= absolute
            hm2 = hm1 / absolute
            hm1 = result
        end
    end

    return result, exponent
end

""" Logarithm of the factorial """
global factorialLog = Vector{Float64}(undef, 2000)
factorialLog[1] = 0
for i in 2:length(factorialLog)
    factorialLog[i] = factorialLog[i - 1] + log(i)
end

function FactorialLog(n)
    global factorialLog

    if n < 0
        return NaN
    elseif n == 0
        return 0.0
    elseif n <= length(factorialLog)
        return factorialLog[n]
    else
        return NaN
    end
end

function WaveFunction(n::Int, x::Real; s=1.0)
    ξ = x * s

    normLog = 0.5 * log(s / sqrt(pi)) - 0.5 * (FactorialLog(n) + n * log(2.0))

    hermite, exponent = HermiteLog(n, ξ)
    
    result = 0.0

    if hermite != 0.0
        resultLog = log(abs(hermite))
        result = normLog - 0.5 * ξ * ξ + resultLog + exponent
        result = sign(hermite) * exp(result)
    end

    return result
end

function WaveFunctionSquare(ρ::Operator, x; s=1.0)
    result = 0.0
    data = real(ρ.data)

    len = size(data)[1]

    wf = [WaveFunction(n - 1, x, s=s) for n = 1:len]

    for m = 1:len
        for n = 1:len
            result += data[m, n] * wf[n] * wf[m]
        end
    end

    return result
end

function WaveFunction(Ψ::Ket, x; s=1.0)
    result = 0.0 + 0.0im
    data = Ψ.data

    len = length(data)

    for n = 1:len
        result += data[n] * WaveFunction(n - 1, x, s=s)
    end

    return result
end