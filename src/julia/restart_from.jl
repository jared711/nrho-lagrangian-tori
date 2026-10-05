# Turn a torus written by bin/param (output_torus<eps>, possibly .gz) into a param input file, to
# restart or extend a continuation from it. Works for any grid shape (e.g. 128 x 512).
#
# Usage:
#   julia --project=. src/julia/restart_from.jl <output_torus file> <new start.csv> <step>
# <step> is the continuation step written into the input (positive: outward, negative: inward).
# epsilon is reset to 0, so epsilon in the new run is relative to this torus.

function main(args)
    src, dst, step = args[1], args[2], args[3]
    L = endswith(src, ".gz") ? readlines(`gzip -dc $src`) : readlines(src)
    h = L[1:12]                    # toltail tolinva tolinte w1 w2 eps H mu n1 n2 deps 1
    h[6] = "0.0"
    h[11] = step
    open(dst, "w") do f
        println(f, join(h, " "))
        for line in L[13:end]
            v = split(line)
            println(f, parse(Int, v[1]) + 1, " ", parse(Int, v[2]) + 1, " ", join(v[3:6], " "))
        end
    end
    println("wrote $dst: grid $(h[9]) x $(h[10]), step $step")
end

main(ARGS)
