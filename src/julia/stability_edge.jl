# Size of the regular (KAM) region around an elliptic NRHO: how far from the NRHO orbits stay
# quasi-periodic before they become chaotic or escape. Used to choose NRHO members with large tori.
#
# Usage:
#   julia --project=. src/julia/stability_edge.jl <Id> [<Id> ...]
# For each member: fixed point of the section map (q3 = 0) and its elliptic eigendirections V1, V2;
# then orbits started at z* + a (c1 Re V1 + c2 Re V2) for the directions (c1, c2) = (1, 1),
# (1, 0.3), (0.3, 1) and growing a. An orbit counts as regular if the NAFF rotation numbers of the
# two halves of 4096 iterates agree to 1e-7. Reported per direction: the largest regular a, the
# maximum distance of that orbit from the NRHO on the section, and the first a that is not regular.

using CSV, DataFrames, LinearAlgebra, Printf, OrdinaryDiffEq, DiffEqBase, ThreeBodyProblem
include(joinpath(@__DIR__, "util.jl"))
include(joinpath(@__DIR__, "dynamics.jl"))
include(joinpath(@__DIR__, "parameterization-method.jl"))
include(joinpath(@__DIR__, "pmap_tools.jl"))

const LU_KM = 238042.0

function edge(id; dirs=((1.0, 1.0), (1.0, 0.3), (0.3, 1.0)), a0=2e-5, factor=1.15, amax=5e-4)
    root = joinpath(@__DIR__, "..", "..")
    df = CSV.read(joinpath(root, "data", "initial_conditions", "NRHO_L2.csv"), DataFrame, normalizenames=true)
    row = df[df.Id .== id, :][1, :]
    mu, H = row.Mass_ratio, -row.Jacobi_constant_LU2_TU2_ / 2
    rv = [row.x0_LU_, row.y0_LU_, row.z0_LU_, row.vx0_LU_TU_, row.vy0_LU_TU_, row.vz0_LU_TU_]
    x, _, _ = P(rv2pq(rv), mu)
    base = tempname() * "_base.csv"
    write_input(base, [0.1, 0.2], repeat(reshape(x[[1, 2, 4, 5]], 1, 1, 4), 2, 2, 1), H, mu)
    z, D = redirect_stdout(devnull) do
        fixed_point(base, x[[1, 2, 4, 5]])
    end
    all(abs.(abs.(eigvals(D)) .- 1) .< 1e-6) || return @printf("Id %d: not elliptic-elliptic\n", id)
    _, V = eigen_frame(D)
    for (c1, c2) in dirs
        a, last_regular, rad_regular = a0, 0.0, 0.0
        while a <= amax
            pts, _ = run_orbit(base, z + a * c1 * real(V[:, 1]) + a * c2 * real(V[:, 2]), 4096)
            regular = size(pts, 1) == 4096 &&
                      maximum(abs.(refined_rotation_numbers(pts[1:2048, :], z, V) -
                                   refined_rotation_numbers(pts[2049:end, :], z, V))) < 1e-7
            regular || break
            last_regular = a
            rad_regular = maximum(sqrt.((pts[:, 1] .- z[1]) .^ 2 .+ (pts[:, 2] .- z[2]) .^ 2)) * LU_KM
            a *= factor
        end
        @printf("Id %3d  direction (%.1f, %.1f): regular up to a = %.2e (%.1f km from the NRHO on the section); first non-regular a = %.2e\n",
                id, c1, c2, last_regular, rad_regular, a)
    end
    rm(base; force=true)
end

for id in parse.(Int, ARGS)
    edge(id)
end
