# Table and figures for a family of tori computed by bin/param (continuation runs).
#
# Usage:
#   julia --project=. src/julia/family_analysis.jl <outprefix> <label1>=<dir1> [<label2>=<dir2> ...]
# Each directory holds the output_torus<eps> files of one continuation run and its run.log (several
# directories of one branch can be given as dir1,dir2)
# (or fam.log). Writes <outprefix>_table.csv and <outprefix>_*.png.

using Printf, LinearAlgebra, FFTW, Plots
gr()

const LU_KM = 238042.0   # Saturn-Enceladus distance (km), the length unit of the CR3BP

# Validated categorical slots (dataviz reference palette, light mode): blue, orange
const SERIES = ["#2a78d6", "#eb6834"]
const INK = "#3d3d3a"

"""readlines that also accepts gzip-compressed files (.gz)."""
readlines_any(path) = endswith(path, ".gz") ? readlines(`gzip -dc $path`) : readlines(path)

"""Read an output_torus file written by param: 12 header lines, then `i1 i2 q1 q2 p1 p2`."""
function read_output_torus(path)
    L = readlines_any(path)
    h = parse.(Float64, L[1:12])
    n1, n2 = Int(h[9]), Int(h[10])
    K = zeros(n1, n2, 4)
    for line in L[13:end]
        x = parse.(Float64, split(line))
        K[Int(x[1])+1, Int(x[2])+1, :] = x[3:6]
    end
    return (omega=h[4:5], eps=h[6], H=h[7], mu=h[8], K=K)
end

"""Error of invariance of each accepted torus, from the lines `eps n1 n2 twist s1 s2 error` of the log."""
function accepted_errors(logfile)
    errs = Dict{Float64,Float64}()
    for line in readlines_any(logfile)
        v = split(line)
        length(v) == 7 || continue
        x = tryparse.(Float64, v)
        any(isnothing, x) && continue
        errs[round(x[1], digits=6)] = x[7]
    end
    return errs
end

function branch_table(dir)
    files = filter(f -> startswith(f, "output_torus"), readdir(dir))
    logf = first(filter(isfile, joinpath.(dir, ["run.log", "run.log.gz", "fam.log"])))
    errs = accepted_errors(logf)
    rows = []
    for f in files
        t = read_output_torus(joinpath(dir, f))
        c = sum(t.K, dims=(1, 2)) / (size(t.K, 1) * size(t.K, 2))
        dq = t.K[:, :, 1:2] .- c[:, :, 1:2]
        radius = maximum(sqrt.(dq[:, :, 1] .^ 2 .+ dq[:, :, 2] .^ 2))   # max distance from the center in (q1, q2)
        push!(rows, (eps=t.eps, omega1=t.omega[1], omega2=t.omega[2], radius_km=radius * LU_KM,
                     grid=size(t.K, 1), error=get(errs, round(t.eps, digits=6), NaN), K=t.K))
    end
    sort!(rows, by=r -> r.omega2)
    return rows
end

function main(args)
    outprefix = args[1]
    branches = [(split(a, "=")[1], split(a, "=")[2]) for a in args[2:end]]
    # a label may combine several directories (consecutive pieces of one branch), comma-separated;
    # epsilon in the table is relative to the start torus of each piece
    tables = [(lab, sort(vcat([branch_table(d) for d in split(dir, ",")]...), by=r -> r.omega2)) for (lab, dir) in branches]

    open(outprefix * "_table.csv", "w") do io
        println(io, "branch,eps,omega1,omega2,radius_km,grid,error_LU")
        for (lab, rows) in tables, r in rows
            @printf(io, "%s,%.8f,%.15f,%.15f,%.4f,%d,%.3e\n", lab, r.eps, r.omega1, r.omega2, r.radius_km, r.grid, r.error)
        end
    end

    common = (framestyle=:box, grid=:y, gridalpha=0.15, foreground_color_axis=INK, legend=:best,
              titlefontsize=11, guidefontsize=10, tickfontsize=9, size=(640, 420), dpi=200)

    # 1. Family in frequency space
    p1 = plot(; xlabel="rotation number ρ₁", ylabel="rotation number ρ₂",
              title="Family of Lagrangian tori around NRHO Id $(get(ENV, "NRHO_ID", "97"))", common...)
    for (i, (lab, rows)) in enumerate(tables)
        plot!(p1, [r.omega1 for r in rows], [r.omega2 for r in rows]; label=lab, color=SERIES[i],
              lw=2, marker=:circle, ms=3, msw=0)
    end
    if haskey(ENV, "RHO_LINEAR")   # linear rotation numbers of the fixed point (the NRHO), "r1,r2"
        rl = parse.(Float64, split(ENV["RHO_LINEAR"], ","))
        scatter!(p1, [rl[1]], [rl[2]]; label="NRHO (linear)", color=INK, marker=:diamond, ms=6, msw=0)
    end
    savefig(p1, outprefix * "_frequencies.png")

    # 2. Torus size along the family
    p2 = plot(; xlabel="rotation number ρ₂", ylabel="torus radius in (q₁, q₂)  [km]",
              title="Torus size along the family", common...)
    for (i, (lab, rows)) in enumerate(tables)
        plot!(p2, [r.omega2 for r in rows], [r.radius_km for r in rows]; label=lab, color=SERIES[i],
              lw=2, marker=:circle, ms=3, msw=0)
    end
    savefig(p2, outprefix * "_size.png")

    # 3. Error of invariance per torus
    p3 = plot(; xlabel="rotation number ρ₂", ylabel="error of invariance  max|E|  [LU]", yscale=:log10,
              title="Accuracy of the computed tori", common...)
    for (i, (lab, rows)) in enumerate(tables)
        ok = [r for r in rows if !isnan(r.error)]
        scatter!(p3, [r.omega2 for r in ok], [r.error for r in ok]; label=lab, color=SERIES[i], ms=4, msw=0)
    end
    savefig(p3, outprefix * "_error.png")

    # 4. Tori in the section plane: small multiples of the smallest, a middle and the largest torus,
    #    each showing every grid point, on a common scale and centered at the smallest torus
    allrows = sort(vcat([rows for (_, rows) in tables]...), by=r -> r.radius_km)
    pick = allrows[unique([1, cld(length(allrows), 2), length(allrows)])]
    c0 = sum(allrows[1].K, dims=(1, 2)) / length(allrows[1].K[:, :, 1])
    lim = 1.1 * maximum(maximum(abs.(r.K[:, :, k] .- c0[k])) for r in pick for k in 1:2) * LU_KM
    panels = [scatter(vec(r.K[:, :, 1] .- c0[1]) .* LU_KM, vec(r.K[:, :, 2] .- c0[2]) .* LU_KM;
                      color=SERIES[1], ms=1.5, msw=0, label="", aspect_ratio=:equal, xlims=(-lim, lim),
                      ylims=(-lim, lim), framestyle=:box, grid=false, title=@sprintf("radius %.1f km", r.radius_km),
                      titlefontsize=10, xlabel="q₁ − q₁* [km]", ylabel=(j == 1 ? "q₂ − q₂* [km]" : ""),
                      guidefontsize=9, tickfontsize=8) for (j, r) in enumerate(pick)]
    p4 = plot(panels...; layout=(1, length(panels)), size=(1000, 380), dpi=200, left_margin=6Plots.mm, bottom_margin=4Plots.mm,
              plot_title="Tori on the section q₃ = 0 (grid points of K)", plot_titlefontsize=11)
    savefig(p4, outprefix * "_section.png")

    for (lab, rows) in tables
        @printf("%-10s %3d tori, rho2 in [%.6f, %.6f], radius %.2f..%.2f km, error median %.2e max %.2e\n",
                lab, length(rows), minimum(r.omega2 for r in rows), maximum(r.omega2 for r in rows),
                minimum(r.radius_km for r in rows), maximum(r.radius_km for r in rows),
                median_([r.error for r in rows if !isnan(r.error)]), maximum([r.error for r in rows if !isnan(r.error)]; init=-Inf))
    end
end

median_(v) = isempty(v) ? NaN : (s = sort(v); n = length(s); isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2)

main(ARGS)
