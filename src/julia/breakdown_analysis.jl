# Limits of the torus families: how far the regular region extends around each NRHO, and how the
# tori break down as they approach its edge.
#
# Usage:
#   julia --project=. src/julia/breakdown_analysis.jl <outprefix>
# Inputs: results/families/stability_edge_scan.txt, elliptic_members_scan.txt, the torus files in
# results/families/, and new orbit sweeps (computed here) for NRHO Id 53, mix (1, 1), and Id 51,
# mix (0.3, 1). Writes <outprefix>_extent.png, _drift.png, _modes.png and CSV files with the data.

using CSV, DataFrames, LinearAlgebra, Printf, FFTW, Plots, OrdinaryDiffEq, DiffEqBase, ThreeBodyProblem
gr()
include(joinpath(@__DIR__, "util.jl"))
include(joinpath(@__DIR__, "dynamics.jl"))
include(joinpath(@__DIR__, "parameterization-method.jl"))
include(joinpath(@__DIR__, "pmap_tools.jl"))

const ROOT = joinpath(@__DIR__, "..", "..")
const FAM = joinpath(ROOT, "results", "families")
const LU_KM = 238042.0
const SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]   # validated categorical slots (light)
const INK = "#3d3d3a"
const COMMON = (framestyle=:box, grid=:y, gridalpha=0.15, titlefontsize=11, guidefontsize=10,
                tickfontsize=9, dpi=200)

readlines_any(path) = endswith(path, ".gz") ? readlines(`gzip -dc $path`) : readlines(path)

# ---- 1. extent of the regular region per member ------------------------------------------------
function extent_table()
    peri = Dict{Int,Float64}()
    for l in readlines(joinpath(FAM, "elliptic_members_scan.txt"))
        v = split(l); length(v) > 5 && v[1] == "Id" && (peri[parse(Int, v[2])] = parse(Float64, v[6]))
    end
    rows = []
    for l in readlines(joinpath(FAM, "stability_edge_scan.txt"))
        m = match(r"Id\s+(\d+)\s+direction \(([\d.]+), ([\d.]+)\): regular up to a = ([\d.e+-]+) \(([\d.]+) km", l)
        m === nothing && continue
        id = parse(Int, m[1])
        push!(rows, (id=id, mix="($(m[2]), $(m[3]))", a=parse(Float64, m[4]), extent=parse(Float64, m[5]),
                     perilune=get(peri, id, NaN)))
    end
    return rows
end

# largest computed torus per member (all families and slices)
function torus_radius(path)
    L = readlines_any(path); h = parse.(Float64, L[1:12]); n = Int(h[9]) * Int(h[10])
    q = zeros(n, 2)
    for (i, line) in enumerate(L[13:end]); v = split(line); q[i, 1] = parse(Float64, v[3]); q[i, 2] = parse(Float64, v[4]); end
    c = sum(q, dims=1) / n
    return maximum(sqrt.(sum((q .- c) .^ 2, dims=2))) * LU_KM
end

function largest_computed()
    best = Dict{Int,Float64}()
    for d in readdir(FAM)
        m = match(r"nrho_id(\d+)", d); m === nothing && continue
        id = parse(Int, m[1])
        for (root, _, files) in walkdir(joinpath(FAM, d)), f in files
            startswith(f, "output_torus") || continue
            best[id] = max(get(best, id, 0.0), torus_radius(joinpath(root, f)))
        end
    end
    return best
end

# ---- 2. regularity sweep across the edge ---------------------------------------------------------
function sweep(id, c1, c2; amps=exp.(range(log(2e-5), log(2.2e-4), length=26)))
    df = CSV.read(joinpath(ROOT, "data", "initial_conditions", "NRHO_L2.csv"), DataFrame, normalizenames=true)
    row = df[df.Id .== id, :][1, :]
    mu, H = row.Mass_ratio, -row.Jacobi_constant_LU2_TU2_ / 2
    x, _, _ = P(rv2pq([row.x0_LU_, row.y0_LU_, row.z0_LU_, row.vx0_LU_TU_, row.vy0_LU_TU_, row.vz0_LU_TU_]), mu)
    base = tempname() * "_base.csv"
    write_input(base, [0.1, 0.2], repeat(reshape(x[[1, 2, 4, 5]], 1, 1, 4), 2, 2, 1), H, mu)
    z, D = redirect_stdout(devnull) do; fixed_point(base, x[[1, 2, 4, 5]]); end
    _, V = eigen_frame(D)
    out = []
    for a in amps
        pts, _ = run_orbit(base, z + a * c1 * real(V[:, 1]) + a * c2 * real(V[:, 2]), 4096)
        if size(pts, 1) < 4096
            push!(out, (a=a, dist=NaN, drift=NaN, escaped=true)); continue
        end
        dist = maximum(sqrt.((pts[:, 1] .- z[1]) .^ 2 .+ (pts[:, 2] .- z[2]) .^ 2)) * LU_KM
        drift = maximum(abs.(refined_rotation_numbers(pts[1:2048, :], z, V) - refined_rotation_numbers(pts[2049:end, :], z, V)))
        push!(out, (a=a, dist=dist, drift=drift, escaped=dist > 1000))   # escaped: leaves Enceladus
    end
    rm(base; force=true)
    return out
end

# ---- 3. Fourier content of the tori along a family -----------------------------------------------
"""Highest |k_1| and |k_2| whose Fourier coefficient (physical units, LU) exceeds `tol`: the modes
needed to represent the torus to the accuracy it is solved to (tol-floor 1e-11 LU)."""
function spectral_width(path; tol=1e-11)
    L = readlines_any(path); h = parse.(Float64, L[1:12]); n1, n2 = Int(h[9]), Int(h[10])
    K = zeros(n1, n2, 4)
    for line in L[13:end]; v = split(line); K[parse(Int, v[1])+1, parse(Int, v[2])+1, :] = parse.(Float64, v[3:6]); end
    A = zeros(n1, n2)
    for c in 1:4
        C = abs.(fft(K[:, :, c] .- sum(K[:, :, c]) / (n1 * n2))) / (n1 * n2)
        A = max.(A, C)
    end
    kk(j, n) = j < n ÷ 2 ? j : j - n
    k1 = 0; k2 = 0
    for i in 1:n1, j in 1:n2
        if A[i, j] > tol
            k1 = max(k1, abs(kk(i - 1, n1))); k2 = max(k2, abs(kk(j - 1, n2)))
        end
    end
    return k1, k2, "$(n1)x$(n2)"
end

function family_modes(dirs)
    rows = []
    for d in dirs, f in readdir(d)
        startswith(f, "output_torus") || continue
        p = joinpath(d, f)
        k1, k2, g = spectral_width(p)
        push!(rows, (radius=torus_radius(p), k1=k1, k2=k2, grid=g))
    end
    sort!(rows, by=r -> r.radius)
end

function main(args)
    pre = args[1]

    # 1.
    ext = extent_table()
    comp = largest_computed()
    CSV.write(pre * "_extent.csv", DataFrame(ext))
    p1 = plot(; xlabel="perilune altitude of the NRHO [km]", ylabel="distance from the NRHO on the section [km]",
              title="Extent of the regular region around each NRHO", legend=:outertopright, size=(820, 440), COMMON...)
    for (i, mix) in enumerate(["(1.0, 1.0)", "(1.0, 0.3)", "(0.3, 1.0)"])
        r = [e for e in ext if e.mix == mix]
        scatter!(p1, [e.perilune for e in r], [e.extent for e in r]; color=SERIES[i], ms=4, msw=0,
                 label="regular extent, mode mix $mix")
    end
    ids = sort(collect(keys(comp)))
    pmap = Dict(e.id => e.perilune for e in ext)
    scatter!(p1, [pmap[i] for i in ids if haskey(pmap, i)], [comp[i] for i in ids if haskey(pmap, i)];
             color=INK, marker=:diamond, ms=6, msw=0, label="largest computed torus (radius)")
    for i in ids
        haskey(pmap, i) && annotate!(p1, pmap[i] + 0.3, comp[i], text("Id $i", 7, INK, :left))
    end
    savefig(p1, pre * "_extent.png")

    # 2.
    s53 = sweep(53, 1.0, 1.0); s51 = sweep(51, 0.3, 1.0)
    CSV.write(pre * "_drift.csv", vcat(DataFrame(s53) |> d -> (d.member .= "Id 53 (1,1)"; d),
                                       DataFrame(s51) |> d -> (d.member .= "Id 51 (0.3,1)"; d)))
    p2 = plot(; xlabel="maximum distance of the orbit from the NRHO on the section [km]",
              ylabel="frequency drift between orbit halves", yscale=:log10, legend=:topleft,
              title="Regular orbits become chaotic or escape at the edge", size=(820, 440), COMMON...)
    for (i, (lab, s)) in enumerate((("Id 53, mix (1, 1)", s53), ("Id 51, mix (0.3, 1)", s51)))
        ok = [r for r in s if !r.escaped]
        plot!(p2, [r.dist for r in ok], [max(r.drift, 1e-13) for r in ok]; color=SERIES[i], lw=2, marker=:circle,
              ms=3, msw=0, label=lab)
        esc = [r for r in s if r.escaped]
        if !isempty(esc)
            last_ok = maximum(r.dist for r in ok)
            vline!(p2, [last_ok]; color=SERIES[i], ls=:dash, lw=1, label="$lab: orbits escape beyond")
        end
    end
    hline!(p2, [1e-7]; color=INK, ls=:dot, lw=1, label="regularity threshold")
    savefig(p2, pre * "_drift.png")

    # 3.
    m53 = family_modes([joinpath(FAM, "nrho_id53", d) for d in readdir(joinpath(FAM, "nrho_id53")) if isdir(joinpath(FAM, "nrho_id53", d))])
    m51 = family_modes([joinpath(FAM, "nrho_id51", "mix_0.3_1", "outward_1")])
    CSV.write(pre * "_modes.csv", vcat(DataFrame(m53) |> d -> (d.member .= "Id 53"; d), DataFrame(m51) |> d -> (d.member .= "Id 51"; d)))
    p3 = plot(; xlabel="torus radius on the section [km]", ylabel="highest significant Fourier mode |kⱼ|",
              title="Fourier modes needed to 1e-11 LU grow toward the edge",
              legend=:topleft, size=(820, 440), COMMON...)
    for (i, (lab, m)) in enumerate((("Id 53", m53), ("Id 51", m51)))
        plot!(p3, [r.radius for r in m], [r.k1 for r in m]; color=SERIES[i], lw=2, label="$lab, angle θ₁")
        plot!(p3, [r.radius for r in m], [r.k2 for r in m]; color=SERIES[i], lw=2, ls=:dash, label="$lab, angle θ₂")
    end
    savefig(p3, pre * "_modes.png")
    println("done: ", pre)
end

main(ARGS)
