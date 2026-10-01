# Set up the continuation of a family of Lagrangian tori around an elliptic-elliptic NRHO of
# data/initial_conditions/NRHO_L2.csv: fixed point of the section map, linear rotation numbers,
# a regular orbit-fitted start torus, and the continuation direction.
#
# Usage:
#   julia --project=. src/julia/family_pipeline.jl <Id> <outdir>
# Writes <outdir>/start.csv (param input, 32x32), <outdir>/setup.txt (rho, rho_lin, domega and
# the param commands for the two branches) and results/families/nrho_id<Id>/nrho_fixed_point.txt.
# The amplitudes tried for the start torus are 4e-5, 2e-5 and 1e-5 LU along both elliptic
# eigendirections; the first whose orbit is regular (NAFF frequency drift between the two halves of
# 8192 iterates below 1e-8) and whose Fourier fit residual is below 1e-9 is used.

using CSV, DataFrames, LinearAlgebra, Printf, OrdinaryDiffEq, DiffEqBase, ThreeBodyProblem
include(joinpath(@__DIR__, "util.jl"))
include(joinpath(@__DIR__, "dynamics.jl"))
include(joinpath(@__DIR__, "parameterization-method.jl"))   # P(): section map of the 6D flow
include(joinpath(@__DIR__, "pmap_tools.jl"))

function main(args)
    id, outdir = parse(Int, args[1]), args[2]
    mkpath(outdir)
    root = joinpath(@__DIR__, "..", "..")
    df = CSV.read(joinpath(root, "data", "initial_conditions", "NRHO_L2.csv"), DataFrame, normalizenames=true)
    row = df[df.Id .== id, :][1, :]
    mu, H, T = row.Mass_ratio, -row.Jacobi_constant_LU2_TU2_ / 2, row.Period_TU_
    rv = [row.x0_LU_, row.y0_LU_, row.z0_LU_, row.vx0_LU_TU_, row.vy0_LU_TU_, row.vz0_LU_TU_]

    # fixed point of the section map (q3 = 0, upward crossing)
    x, _, _ = P(rv2pq(rv), mu)
    base = joinpath(outdir, "base.csv")
    write_input(base, [0.1, 0.2], repeat(reshape(x[[1, 2, 4, 5]], 1, 1, 4), 2, 2, 1), H, mu)
    z, D = fixed_point(base, x[[1, 2, 4, 5]])
    ev = eigvals(D)
    all(abs.(abs.(ev) .- 1) .< 1e-6) || error("Id $id is not elliptic-elliptic: |lambda| = $(abs.(ev))")
    rho_lin, V = eigen_frame(D)
    famdir = joinpath(root, "results", "families", "nrho_id$id")
    mkpath(famdir)
    open(joinpath(famdir, "nrho_fixed_point.txt"), "w") do f
        println(f, join(string.(vcat(z, H, mu, T)), " "))
    end

    # regular orbit-fitted start torus
    chosen = nothing
    for a in (4e-5, 2e-5, 1e-5)
        pts, _ = run_orbit(base, z + a * real(V[:, 1]) + a * real(V[:, 2]), 8192)
        size(pts, 1) < 8192 && continue
        r = refined_rotation_numbers(pts, z, V)
        drift = maximum(abs.(refined_rotation_numbers(pts[1:4096, :], z, V) - refined_rotation_numbers(pts[4097:end, :], z, V)))
        K, res = fit_torus(pts[1:6000, :], r, 10, 32)
        @printf("Id %d  a = %.0e: rho = (%.9f, %.9f), drift %.1e, fit residual %.1e\n", id, a, r[1], r[2], drift, res)
        if drift < 1e-8 && res < 1e-9
            chosen = (a=a, rho=r, K=K)
            break
        end
    end
    chosen === nothing && error("no regular start torus found for Id $id")
    domega = chosen.rho - rho_lin
    write_input(joinpath(outdir, "start.csv"), chosen.rho, chosen.K, H, mu)
    open(joinpath(outdir, "setup.txt"), "w") do f
        @printf(f, "Id %d  H = %.15g  mu = %.15g  T = %.15g\n", id, H, mu, T)
        @printf(f, "rho_lin = %.17g %.17g\nrho = %.17g %.17g\ndomega = %.6e %.6e\namplitude = %.0e\n",
                rho_lin..., chosen.rho..., domega..., chosen.a)
        @printf(f, "inward:  param start.csv --domega %.6e %.6e --eps-max -0.95 --max-steps 3000 --tol-floor 1e-11  (step -0.01)\n", domega...)
        @printf(f, "outward: param start.csv --domega %.6e %.6e --eps-max 100 --max-steps 3000 --tol-floor 1e-11  (step +0.01)\n", domega...)
    end
    print(read(joinpath(outdir, "setup.txt"), String))
end

main(ARGS)
