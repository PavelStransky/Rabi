""" Matrix export (for Origin) """
function Export(fname::String, xs, ys=nothing)
    fn = fname * ".txt"
    open(fn, "w") do f
        for i = 1:length(xs)
            write(f, "$(xs[i])\t")
            if ys !== nothing 
                write(f, join(ys[i,:], "\t")) 
            end
            write(f, "\n")
        end
    end
    println("Saved into $fn")
end

" Matrix export (for Origin and GCM) "
function Export(fname::String, xs, ys, data::Array{Float64,3})
    num = size(data)[1]

    for k = 1:num
        fn = "$(fname)_$k.txt"
        open(fn, "w") do f
            for i = 1:length(xs)
                for j = 1:length(ys)
                    write(f, "$(xs[i])\t$(ys[j])\t$(data[k,j,i])\n")
                end
            end
        end
        println("Saved into $fn")

        fn = "$(fname)_gcm_$k.txt"
        open(fn, "w") do g
            write(g, "PavelStransky.Math.Matrix\n")
            write(g, "$(length(xs))\t$(length(ys))\n")
            for i = 1:length(xs)
                write(g, join(data[k,:,i], "\t"))
                write(g, "\n")
            end
        end
        println("Saved into $fn")
    end
end

" Export of the matrix (for Origin) "
function Export(fname::String, xs, ys, zs::Array{Float64,2})
    open(fname, "w") do f
        for i = 1:length(xs)
            for j = 1:length(ys)
                write(f, "$(xs[i])\t$(ys[j])\t$(zs[i,j])\n")
            end
        end
    end
end

" Export of the time series (Origin) "
function Export(fname, ts, xs, xm)
    open(fname, "w") do f
        if !isnothing(xm)
            write(f, "-1\t$xm\n")
        end
        for (t, x) in zip(ts, xs)
            write(f, "$t\t$x\n")
        end
    end
end