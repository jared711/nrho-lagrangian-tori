# Helpers that drive the C++ Poincare map (bin/param --orbit) from Julia.
#
# Used to (1) Newton-solve for the fixed point of the section map (the NRHO),
# (2) iterate orbits for frequency analysis, and (3) write initial tori for param.

using LinearAlgebra
using Printf

const PARAM = joinpath(@__DIR__, "..", "..", "bin", "param")

"""
    run_orbit(input_file, z, niter) -> (pts, D)

Iterate the Poincare map `niter` times from `z = (q1,q2,p1,p2)`.
Returns the iterates (niter x 4) and the Jacobian of the map at `z`.
"""
function run_orbit(input_file, z, niter)
    args = [@sprintf("%.17g", v) for v in z]
    out = readlines(`$PARAM $input_file --orbit $niter $args`)
    D = nothing
    pts = Vector{Vector{Float64}}()
    for line in out
        if startswith(line, "# Dfz")
            D = permutedims(reshape(parse.(Float64, split(line)[3:end]), 4, 4))
        elseif startswith(line, "# map failed")
            break
        elseif !startswith(line, "#") && length(split(line)) == 4
            v = tryparse.(Float64, split(line))
            all(!isnothing, v) && push!(pts, v)
        end
    end
    return reduce(vcat, permutedims.(pts)), D
end

"""
    fixed_point(input_file, z0; tol=1e-13, maxit=20) -> (z, D)

Newton's method for P(z) = z.
"""
function fixed_point(input_file, z0; tol=1e-13, maxit=20)
    z = collect(Float64, z0)
    D = nothing
    for it in 0:maxit-1
        pts, D = run_orbit(input_file, z, 1)
        r = pts[1, :] - z
        @printf("  fixed-point it %d |P(z)-z| = %.3e\n", it, norm(r))
        norm(r) < tol && break
        z -= (D - I) \ r
    end
    return z, D
end

"""
    eigen_frame(D) -> (rho, V)

Rotation numbers in (0, 1/2) and complex eigenvectors (as columns), one per
conjugate pair, of the Jacobian at an elliptic-elliptic fixed point.
"""
function eigen_frame(D)
    F = eigen(D)
    keep = findall(λ -> imag(λ) > 0, F.values)
    rho = [angle(F.values[i]) / 2π for i in keep]
    return rho, F.vectors[:, keep]
end

"""
    rotation_numbers(pts, zfix, V) -> (rho, |w|)

Rotation numbers of an orbit, from a linear fit of the unwrapped arguments of its
complex eigen-coordinates w = V^{-1} (z - zfix) against the iterate number.
"""
function rotation_numbers(pts, zfix, V)
    w = ([V conj(V)] \ permutedims(pts .- zfix'))[1:2, :]
    k = collect(0:size(pts, 1)-1)
    rho = zeros(2)
    for j in 1:2
        ang = unwrap_angle(angle.(w[j, :]))
        A = [k ones(length(k))]
        rho[j] = (A \ ang)[1] / 2π
    end
    return rho, abs.(w)
end

function unwrap_angle(a)
    out = copy(a)
    for i in 2:length(a)
        d = a[i] - a[i-1]
        out[i] = out[i-1] + (d - 2π * round(d / 2π))
    end
    return out
end

"""
    linear_torus(zfix, V, amps, n) -> K (n x n x 4)

First-order torus K(θ) = zfix + Σ_j a_j Re(V_j exp(2πi θ_j)) on an n x n grid.
"""
function linear_torus(zfix, V, amps, n)
    K = zeros(n, n, 4)
    for i1 in 0:n-1, i2 in 0:n-1
        e1, e2 = cis(2π * i1 / n), cis(2π * i2 / n)
        K[i1+1, i2+1, :] = zfix + amps[1] * real(V[:, 1] * e1) + amps[2] * real(V[:, 2] * e2)
    end
    return K
end

"""
    write_input(path, rho, K, H, mu; ...)

Write a param input file. `K[i1, i2, :] = K(i1/n, i2/n)`; the second index runs
fastest, matching position()/indices() in grid.h.
"""
function write_input(path, rho, K, H, mu; toltail=1e-11, tolinva=1e-12, tolinte=1e-17, deps=5e-5)
    n = size(K, 1)
    open(path, "w") do f
        @printf(f, "%.17g %.17g %.17g %.17g %.17g 0.0 %.17g %.17g %d %d %.17g 1\n",
                toltail, tolinva, tolinte, rho[1], rho[2], H, mu, n, n, deps)
        for i1 in 1:n, i2 in 1:n
            println(f, i1, " ", i2, " ", join([@sprintf("%.17g", v) for v in K[i1, i2, :]], " "))
        end
    end
end

"""
    naff_frequency(w, guess; width) -> ν

Refine the dominant frequency of the complex signal w[k] near `guess` (cycles per
iterate) by maximizing the Hann-windowed Fourier amplitude (Laskar's NAFF).
"""
function naff_frequency(w, guess; width=2 / length(w))
    N = length(w)
    k = 0:N-1
    χ = 1 .- cos.(2π .* k ./ N)
    amp(ν) = -abs(sum(χ .* w .* cis.(-2π * ν .* k)))
    a, b = guess - width, guess + width
    φ = (sqrt(5) - 1) / 2
    c, d = b - φ * (b - a), a + φ * (b - a)
    fc, fd = amp(c), amp(d)
    for _ in 1:80
        if fc < fd
            b, d, fd = d, c, fc
            c = b - φ * (b - a); fc = amp(c)
        else
            a, c, fc = c, d, fd
            d = a + φ * (b - a); fd = amp(d)
        end
    end
    return (a + b) / 2
end

"""
    refined_rotation_numbers(pts, zfix, V) -> rho

Linear-fit rotation numbers refined by `naff_frequency` on each eigen-coordinate.
"""
function refined_rotation_numbers(pts, zfix, V)
    rho0, _ = rotation_numbers(pts, zfix, V)
    w = ([V conj(V)] \ permutedims(pts .- zfix'))[1:2, :]
    return [naff_frequency(w[j, :], rho0[j]) for j in 1:2]
end

"""
    fit_torus(pts, rho, M, n) -> (K, residual)

Least-squares fit of a Fourier torus to an orbit, z_k ≈ K(k ρ) with
K(θ) = Σ_{|m_i| ≤ M} c_m exp(2πi m·θ), then evaluation of K on an n x n grid.
The orbit point z_0 sits at θ = 0.
"""
function fit_torus(pts, rho, M, n)
    N = size(pts, 1)
    modes = [(m1, m2) for m1 in -M:M for m2 in -M:M]
    A = [cis(2π * (m[1] * rho[1] + m[2] * rho[2]) * k) for k in 0:N-1, m in modes]
    C = A \ complex.(pts)
    residual = maximum(abs.(real.(A * C) .- pts))
    K = zeros(n, n, 4)
    for i1 in 0:n-1, i2 in 0:n-1
        e = [cis(2π * (m[1] * i1 + m[2] * i2) / n) for m in modes]
        K[i1+1, i2+1, :] = real.(transpose(e) * C)
    end
    return K, residual
end

"""
    map_noise(input_file, z, v; hs) -> Vector of (h, defect)

Smoothness test of the section map along direction v: the second-order defect
|P(z+hv) - 2P(z) + P(z-hv)| / h^2 should be constant (≈ |D²P(v,v)|) for a smooth
map. It grows like noise/h^2 once h is small enough that integration noise dominates.
"""
function map_noise(input_file, z, v; hs=10.0 .^ (-3:-0.5:-9))
    P0 = run_orbit(input_file, z, 1)[1][1, :]
    out = Tuple{Float64,Float64}[]
    for h in hs
        Pp = run_orbit(input_file, z + h * v, 1)[1][1, :]
        Pm = run_orbit(input_file, z - h * v, 1)[1][1, :]
        push!(out, (h, norm(Pp - 2P0 + Pm)))
    end
    return out
end

"""
    frequency_map(input_file, zfix, V, amps1, amps2; niter=2048)

Laskar-style frequency map: for each pair of amplitudes (a1, a2) start at
zfix + a1 Re V1 + a2 Re V2, iterate the map, and measure the rotation numbers with
`refined_rotation_numbers` on the whole orbit. Also returns a regularity indicator:
the max change of the rotation numbers between the first and second halves of the
orbit (small for orbits on KAM tori, large for chaotic or escaping ones; Inf if the
map failed).
"""
function frequency_map(input_file, zfix, V, amps1, amps2; niter=2048)
    rho = fill(NaN, length(amps1), length(amps2), 2)
    diff = fill(Inf, length(amps1), length(amps2))
    for (i, a1) in enumerate(amps1), (j, a2) in enumerate(amps2)
        pts, _ = run_orbit(input_file, zfix + a1 * real(V[:, 1]) + a2 * real(V[:, 2]), niter)
        size(pts, 1) < niter && continue
        h = niter ÷ 2
        r = refined_rotation_numbers(pts, zfix, V)
        r1 = refined_rotation_numbers(pts[1:h, :], zfix, V)
        r2 = refined_rotation_numbers(pts[h+1:end, :], zfix, V)
        rho[i, j, :] = r
        diff[i, j] = maximum(abs.(r1 - r2))
    end
    return rho, diff
end

"""
    eval_torus(input_file) -> (F, DF)

Evaluate the section map and its Jacobian at every grid point of the torus in
`input_file` (bin/param --eval). F is n x n x 4 and DF is n x n x 4 x 4, indexed like
`write_input` (K[i1, i2, :]).
"""
function eval_torus(input_file, n)
    F = zeros(n, n, 4)
    DF = zeros(n, n, 4, 4)
    cnt = 0
    for line in readlines(`$PARAM $input_file --eval`)
        v = split(line)
        length(v) == 21 || continue
        l = tryparse(Int, v[1]); l === nothing && continue
        i1, i2 = divrem(l, n)
        x = parse.(Float64, v[2:end])
        F[i1+1, i2+1, :] = x[1:4]
        DF[i1+1, i2+1, :, :] = permutedims(reshape(x[5:20], 4, 4))
        cnt += 1
    end
    cnt == n^2 || error("map evaluation failed ($cnt of $(n^2) points)")
    return F, DF
end

"""
    shift_matrix(n, w) -> S

Real n x n matrix of the Fourier shift u(θ) -> u(θ + w) on a grid of n points
(trigonometric interpolation; the Nyquist mode is shifted by its real part).
"""
function shift_matrix(n, w)
    k = [j < n ÷ 2 ? j : j - n for j in 0:n-1]
    mult = [abs(kk) == n ÷ 2 ? complex(cos(2π * kk * w)) : cis(2π * kk * w) for kk in k]
    Fm = [cis(-2π * j * kk / n) for kk in 0:n-1, j in 0:n-1]
    return real(Fm' * Diagonal(mult) * Fm / n)
end

"""
    collocation_newton(K, rho, H, mu; iters=6, tmpfile) -> (K, errors)

Reference solver: full Newton for F(K(θ)) = K(θ + ρ) on the n x n grid with fixed
frequency, solving the dense linearized equation
    DF(K(θ_p)) ΔK(θ_p) - ΔK(θ_p + ρ) = -E(θ_p)
together with the phase conditions Σ_p ∂_{θ_j}K(θ_p) · ΔK(θ_p) = 0 (least squares).
"""
function collocation_newton(K, rho, H, mu; iters=6, tmpfile="/tmp/colloc_torus.csv")
    n = size(K, 1)
    S1, S2 = shift_matrix(n, rho[1]), shift_matrix(n, rho[2])
    S = kron(S1, S2)                 # acts on vectors with the second index fastest
    D1 = kron(shift_derivative(n), Matrix(I, n, n))
    D2 = kron(Matrix(I, n, n), shift_derivative(n))
    flat(A) = vec(permutedims(A, (2, 1)))     # i2 fastest
    unflat(v) = permutedims(reshape(v, n, n), (2, 1))
    errs = Float64[]
    for it in 1:iters
        write_input(tmpfile, rho, K, H, mu)
        F, DF = eval_torus(tmpfile, n)
        E = zeros(n * n, 4)
        for c in 1:4
            E[:, c] = flat(F[:, :, c]) - S * flat(K[:, :, c])
        end
        push!(errs, maximum(abs.(E)))
        @printf("  collocation Newton it %d  max|E| = %.3e\n", it, errs[end])
        it == iters && break
        m = n * n
        A = zeros(4m + 2, 4m)
        for c in 1:4, d in 1:4
            A[(c-1)*m+1:c*m, (d-1)*m+1:d*m] = Diagonal(flat(DF[:, :, c, d]))
        end
        for c in 1:4
            A[(c-1)*m+1:c*m, (c-1)*m+1:c*m] -= S
            A[4m+1, (c-1)*m+1:c*m] = (D1 * flat(K[:, :, c]))'
            A[4m+2, (c-1)*m+1:c*m] = (D2 * flat(K[:, :, c]))'
        end
        b = [-vec(E); 0.0; 0.0]
        dK = A \ b
        for c in 1:4
            K[:, :, c] += unflat(dK[(c-1)*m+1:c*m])
        end
    end
    return K, errs
end

"""Spectral derivative matrix d/dθ on n grid points (θ ∈ [0,1))."""
function shift_derivative(n)
    k = [j < n ÷ 2 ? j : (j == n ÷ 2 ? 0 : j - n) for j in 0:n-1]
    Fm = [cis(-2π * j * kk / n) for kk in 0:n-1, j in 0:n-1]
    return real(Fm' * Diagonal(2π * im .* k) * Fm / n)
end
