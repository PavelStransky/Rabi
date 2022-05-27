" Data import "
function ImportSection(fname)
    data = Array{Float64}(undef, 0)
    numLines = 0

    open(fname, "r") do f
        for lines in readlines(f)
            append!(data, map(x->parse(Float64, x), split(lines, "\t")))
            numLines += 1
        end
    end

    return reshape(data, :, numLines)
end

" Data export "
function ExportSection(fname, λs, averages)
    open(fname, "w") do f
        for i = 1:length(λs)
            write(f, "$(λs[i])\t")
            first = true
            for j = 1:length(averages)
                if first
                    first = false
                else
                    write(f, "\t")
                end
                write(f, string(averages[j][i]))
            end
            write(f, "\n")
        end
    end
end

function MovingAverage(data, width)
    result = copy(data)

    len = length(data)
    shift = -(width ÷ 2) - 1

    for i = 1:len
        num = 0.0
        sum = 0.0
        for j = 1:width
            k = i + j - shift
            if k >= 1 && k <=len
                sum += data[k]
                num += 1
            end
        end
        result[i] = sum / num
    end

    return result
end

function Smooth(data, num)
    result = copy(data)
    len = length(data)

    if len <= 1
        return result
    end
       
    for k = 1:num
        result[1] = (3.0 * data[1] + data[2]) / 4.0
        for i = 2:(len - 1)
            result[i] = (data[i - 1] + 2.0 * data[i] + data[i + 1]) / 4.0
            result[len] = (data[len - 1] + 3.0 * data[len]) / 4.0
        end
        data = copy(result)
    end

    return result
end

data = ImportSection("d:/results/rabi/tmp/mapL_(100,0.5,0.4,0).txt")

x = data[1,2:end]
y = [Smooth(data[i,2:end], 2000) for i = 2:8]

ExportSection("d:/results/rabi/mapL_(100,0.5,0.4,0)S.txt", x, y)

p1 = plot(x, data[5,:])
p1 = plot!(p1, x, y[4,:])
display(plot(p1))