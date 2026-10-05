# Quasi-periodic orbits (QPOs) on computed tori: integrate the full 3-DOF flow from points of
# selected tori and plot them in Enceladus-centred coordinates of the standard rotating frame.
#
# Usage:
#   julia --project=. src/julia/qpo_visualization.jl <outprefix> <nrevs> <torus1> [<torus2> ...]
# The tori are output_torus files written by bin/param. The NRHO (the fixed point of the section
# map) is read from results/families/nrho_id<NRHO_ID>/nrho_fixed_point.txt (env NRHO_ID, default 97; or NRHO_FIXED).

using LinearAlgebra, Printf, OrdinaryDiffEq, Plots, ThreeBodyProblem
gr()
include(joinpath(@__DIR__, "util.jl"))       # pq2rv (Barcelona -> standard rotating frame)
include(joinpath(@__DIR__, "dynamics.jl"))   # CR3BPdynamicsBar

const LU_KM = 238042.0          # Saturn-Enceladus distance
const R_ENCELADUS_KM = 252.1
const SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]   # validated categorical slots (light)
const INK = "#3d3d3a"

"""readlines that also accepts gzip-compressed files (.gz)."""
readlines_any(path) = endswith(path, ".gz") ? readlines(`gzip -dc $path`) : readlines(path)

function read_output_torus(path)
    L = readlines_any(path)
    h = parse.(Float64, L[1:12])
    n1, n2 = Int(h[9]), Int(h[10])
    K = zeros(n1, n2, 4)
    for line in L[13:end]
        x = parse.(Float64, split(line))
        K[Int(x[1])+1, Int(x[2])+1, :] = x[3:6]
    end
    return (omega=h[4:5], H=h[7], mu=h[8], K=K)
end

"""6D Barcelona state on the section q3 = 0 from z = (q1, q2, p1, p2), with p3 > 0 from H."""
function section_to_state(z, H, mu)
    q1, q2, p1, p2 = z
    r1 = sqrt((q1 - mu)^2 + q2^2)
    r2 = sqrt((q1 - mu + 1)^2 + q2^2)
    p3sq = 2 * (H - q2 * p1 + q1 * p2 + (1 - mu) / r1 + mu / r2) - p1^2 - p2^2
    return [q1, q2, 0.0, p1, p2, sqrt(p3sq)]
end

"""Integrate the flow and return Enceladus-centred positions (km, standard rotating frame)."""
function trajectory_km(x0, mu, tf; npts=20000)
    f!(dx, x, p, t) = (dx .= CR3BPdynamicsBar(x, p, t))
    sol = solve(ODEProblem(f!, x0, (0.0, tf), mu), Vern9(); abstol=1e-13, reltol=1e-13,
                saveat=range(0, tf, length=npts))
    R = zeros(length(sol.t), 3)
    for (k, x) in enumerate(sol.u)
        rv = pq2rv(x)
        R[k, :] = (rv[1:3] .- [1 - mu, 0, 0]) .* LU_KM
    end
    return sol.t, R
end

"""Crossings of the plane y = 0 with z > zmin (near apolune), linearly interpolated."""
function apolune_crossings(R, zmin)
    pts = Vector{Vector{Float64}}()
    for k in 1:size(R, 1)-1
        if R[k, 2] * R[k+1, 2] < 0 && R[k, 3] > zmin
            s = R[k, 2] / (R[k, 2] - R[k+1, 2])
            push!(pts, R[k, :] .+ s .* (R[k+1, :] .- R[k, :]))
        end
    end
    return reduce(vcat, permutedims.(pts))
end

function main(args)
    outprefix, nrevs = args[1], parse(Float64, args[2])
    files = args[3:end]
    nrho_id = get(ENV, "NRHO_ID", "97")
    fpfile = get(ENV, "NRHO_FIXED", joinpath(@__DIR__, "..", "..", "results", "families", "nrho_id" * nrho_id, "nrho_fixed_point.txt"))
    fp = vec(parse.(Float64, split(read(fpfile, String))))
    zfix, H, mu, Tnrho = fp[1:4], fp[5], fp[6], fp[7]
    tn, Rn = trajectory_km(section_to_state(zfix, H, mu), mu, Tnrho; npts=40000)

    qpos = []
    for f in files
        t = read_output_torus(f)
        K = t.K
        c = sum(K, dims=(1, 2)) / (size(K, 1) * size(K, 2))
        radius = maximum(sqrt.((K[:, :, 1] .- c[1]) .^ 2 .+ (K[:, :, 2] .- c[2]) .^ 2)) * LU_KM
        tq, Rq = trajectory_km(section_to_state(K[1, 1, :], t.H, mu), mu, nrevs * Tnrho; npts=round(Int, 4000 * nrevs))
        # distance to the NRHO: nearest point of the (densely sampled) periodic orbit
        d = [minimum(sum((Rn .- Rq[k:k, :]) .^ 2, dims=2))^0.5 for k in 1:size(Rq, 1)]
        push!(qpos, (file=f, radius=radius, t=tq ./ Tnrho, R=Rq, d=d))
        alt = minimum(sqrt.(sum(Rq .^ 2, dims=2))) - R_ENCELADUS_KM   # positions are Enceladus-centred
        @printf("%s: section radius %.2f km, distance to NRHO %.2f .. %.2f km, minimum altitude %+.2f km over %.0f revolutions\n",
                basename(f), radius, minimum(d), maximum(d), alt, nrevs)
    end
    sort!(qpos, by=q -> q.radius)

    common = (framestyle=:box, titlefontsize=11, guidefontsize=10, tickfontsize=9, dpi=200)

    # 1. 3D view: NRHO, the largest QPO, Enceladus
    big = qpos[end]
    p1 = plot(Rn[:, 1], Rn[:, 2], Rn[:, 3]; color=INK, lw=1.5, label="NRHO (Id $nrho_id)",
              xlabel="x [km]", ylabel="y [km]", zlabel="z [km]", size=(700, 620),
              title=@sprintf("QPO on a torus of radius %.1f km (%d revolutions)", big.radius, round(Int, nrevs)), common...)
    plot!(p1, big.R[:, 1], big.R[:, 2], big.R[:, 3]; color=SERIES[1], lw=0.4, alpha=0.6, label="QPO")
    # Enceladus as a light wireframe (GR has no parametric surfaces)
    a = range(0, 2π, length=60)
    for lat in range(-75, 75, step=25) .* (π / 180)
        plot!(p1, R_ENCELADUS_KM .* cos(lat) .* cos.(a), R_ENCELADUS_KM .* cos(lat) .* sin.(a), fill(R_ENCELADUS_KM * sin(lat), length(a));
              color="#b5b4ad", lw=0.6, label="")
    end
    for lon in range(0, π, length=7)[1:end-1]
        plot!(p1, R_ENCELADUS_KM .* cos.(a) .* cos(lon), R_ENCELADUS_KM .* cos.(a) .* sin(lon), R_ENCELADUS_KM .* sin.(a);
              color="#b5b4ad", lw=0.6, label="")
    end
    # equal scales on the three axes (a cube around the NRHO), so Enceladus is a sphere
    ctr = [(maximum(Rn[:, k]) + minimum(Rn[:, k])) / 2 for k in 1:3]
    half = 0.55 * maximum(maximum(Rn[:, k]) - minimum(Rn[:, k]) for k in 1:3)
    plot!(p1; xlims=(ctr[1] - half, ctr[1] + half), ylims=(ctr[2] - half, ctr[2] + half), zlims=(ctr[3] - half, ctr[3] + half))
    savefig(p1, outprefix * "_3d.png")

    # 2. Distance to the NRHO along each QPO
    p2 = plot(; xlabel="time [NRHO periods]", ylabel="distance to the NRHO [km]", legend=:outertopright,
              title="QPOs stay in a thin tube around the NRHO", size=(760, 400), grid=:y, gridalpha=0.15, common...)
    for (i, q) in enumerate(qpos)
        lab = @sprintf("torus radius %.1f km", q.radius)
        plot!(p2, q.t, q.d; color=SERIES[i], lw=1, label=lab)
        annotate!(p2, q.t[end], maximum(q.d), text(@sprintf("%.1f km", q.radius), 8, INK, :right, :bottom))
    end
    savefig(p2, outprefix * "_distance.png")

    # 3. Cross-sections of the QPO tubes at apolune (y = 0, z > 0 high), relative to the NRHO crossing
    zmin = 0.5 * maximum(Rn[:, 3])
    cn = apolune_crossings(Rn, zmin)[1, :]
    panels = []
    lim = 0.0
    # longer integrations (ncross revolutions) for the cross-sections: one crossing per revolution
    ncross = parse(Int, get(ENV, "NCROSS", "300"))
    cuts = [apolune_crossings(trajectory_km(section_to_state(read_output_torus(q.file).K[1, 1, :], H, mu), mu,
                                            ncross * Tnrho; npts=2000 * ncross)[2], zmin) for q in qpos]
    for c in cuts
        lim = max(lim, maximum(abs.(c[:, 1] .- cn[1])), maximum(abs.(c[:, 3] .- cn[3])))
    end
    lim *= 1.1
    for (i, (q, c)) in enumerate(zip(qpos, cuts))
        push!(panels, scatter(c[:, 1] .- cn[1], c[:, 3] .- cn[3]; color=SERIES[1], ms=2, msw=0, label="",
                              aspect_ratio=:equal, xlims=(-lim, lim), ylims=(-lim, lim), framestyle=:box, grid=false,
                              title=@sprintf("torus radius %.1f km (%d crossings)", q.radius, size(c, 1)), titlefontsize=10,
                              xlabel="Δx [km]", ylabel=(i == 1 ? "Δz [km]" : ""), guidefontsize=9, tickfontsize=8))
    end
    p3 = plot(panels...; layout=(1, length(panels)), size=(340 * length(panels), 380), dpi=200,
              left_margin=6Plots.mm, bottom_margin=4Plots.mm,
              plot_title="QPO crossings of y = 0 near apolune, relative to the NRHO", plot_titlefontsize=11)
    savefig(p3, outprefix * "_apolune.png")
end

main(ARGS)
