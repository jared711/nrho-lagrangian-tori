using LinearAlgebra
using FFTW
using OrdinaryDiffEq

# Import terminate! - it's in DiffEqBase which is loaded by OrdinaryDiffEq
# but not automatically exported in newer versions
using DiffEqBase: terminate!

# Conventions
# x = [q₁, q₂, q₃, p₁, p₂, p₃]
# z = [q₁, q₂, p₁, p₂]
"""
    p3(z, H, μ)

Returns the third component of the momentum vector given the reduced state vector z, the Hamiltonian H, and the mass parameter μ.
"""
function p3(z, H, μ)
    # H = params.H
    # μ = params.μ
    q₁, q₂, p₁, p₂ = z
    q₃ = 0.0  # We're on the Poincare section q₃ = 0
    r₁ = sqrt((q₁-μ)^2 + q₂^2 + q₃^2)
    r₂ = sqrt((q₁+1-μ)^2 + q₂^2 + q₃^2)
    p₃ = √(2(H - q₂*p₁ + q₁*p₂ + (1-μ)/r₁ + μ/r₂) - p₁^2 - p₂^2)
    if q₂ < 0 # sign convention for p₃, since the √ operator is always positive
        if q₁ < (μ-1) # L₂
            p₃ = -p₃
        end
    elseif q₂ > 0
        if q₁ > (μ-1) # L₁
            p₃ = -p₃
        end
    end
    return p₃
end

"""
    Dp3(z, μ)

Returns the differential of the third component of the momentum vector with respect to the reduced state vector z.
"""
function Dp3(z, μ)
    q₁, q₂, p₁, p₂ = z
    r₁³= ((q₁-μ)^2   + q₂^2)^1.5 # distance to m1, LARGER MASS
    r₂³= ((q₁+1-μ)^2 + q₂^2)^1.5 # distance to m2, smaller mass

    ∂p₃_∂q₁ = 2( p₂ - (1-μ)*(q₁-μ)/r₁³ - μ*(q₁+1-μ)/r₂³)
    ∂p₃_∂q₂ = 2(-p₁ - (1-μ)*    q₂/r₁³ - μ*      q₂/r₂³)
    ∂p₃_∂p₁ = -2q₂ - 2p₁
    ∂p₃_∂p₂ =  2q₁ - 2p₂
    Dp₃ = [∂p₃_∂q₁, ∂p₃_∂q₂, ∂p₃_∂p₁, ∂p₃_∂p₂]

    if q₂ < 0 # sign convention for p₃, since the √ operator is always positive
        if q₁ < (μ-1) # L₂
            Dp₃ = -Dp₃
        end
    elseif q₂ > 0
        if q₁ > (μ-1) # L₁
            Dp₃ = -Dp₃
        end
    end

    return Dp₃
end

"""
    γ(z, H, μ)

Changes the reduced state vector z into the full state vector x.
"""
function γ(z, H, μ) # changes z∈ℝ⁴  into x∈ℝ⁶
    x = zeros(6)
    x[[1,2,4,5]] = z # note q₃ = 0
    x[6] = p3(z, H, μ) # update x[6] = p₃ from the Hamiltonian
    return x
end

"""
    Dγ(z, μ)

Returns the differential of the full state vector x with respect to the reduced state vector z.
"""
Dγ(z, μ) = [1 0 0 0;
            0 1 0 0;
            0 0 0 0;
            0 0 1 0;
            0 0 0 1;
            Dp3(z, μ)']

"""
    γ⁻¹(x)

Changes the full state vector x (which should be on the poincare section Σ) into the reduced state vector z.
"""
function γ⁻¹(x) # changes x ∈ Σ ⊂ ℝ⁶ into z ∈ ℝ⁴
    @assert isapprox(x[3], 0, atol=1e-10) "x[3] must be zero, got $(x[3])"
    # @assert isapprox(computeH(x, μ) ≈ H) "The Hamiltonian is not conserved."
    z = x[[1,2,4,5]]
    return z
end

"""
    Dγ⁻¹(x)

Returns the differential of the reduced state vector z with respect to the full state vector x.
"""
Dγ⁻¹(x) = [1 0 0 0 0 0;
           0 1 0 0 0 0;
           0 0 0 1 0 0;
           0 0 0 0 1 0]


# The symplectic form
Ω(z) = [0  0 -1  0;
        0  0  0 -1;
        1  0  0  0;
        0  1  0  0]
# ω(ξ,η,z) = ξ'*Ω(z)*η

# The metric
G(z) = I + Dp3(z)'*Dp3(z)
# g(ξ,η,z) = ξ'*G(z)*η

# The poincaré map function and differential
σ(x) = x[3]
Dσ = [0,0,1,0,0,0]'


# write documentation
"""
    P(x₀, μ, tmax=10)

Returns the poincare map, the time, and the differential of the poincare map
"""
function P(x₀, μ, tmax=10) # Poincare Map

    # set up ODE probl
    Φ₀ = I(6) # Initialization of the STM, Φ₀ = I
    w₀ = [x₀; reshape(Φ₀,36,1)] # Reshape the matrix into a vector and append it to the state vector
    tspan = (0., tmax) # integrate from 0 to T₀
    prob = ODEProblem(CR3BPstmBar!,w₀,tspan,μ) # CR3BPstm! is our in-place dynamics function for state and STM

    # event function
    σ(x) = x[3]
    Dσ = [0,0,1,0,0,0]' # defining it as an adjoint rather than a 1x6 matrix makes the linear algebra work out better
    condition(u, t, integrator) = σ(u) # event when x = 0
    function affect!(integrator)
        integrator.u[3] = 0.0 # actually set z = 0, to prevent the event from triggering again
        terminate!(integrator)
    end
    cb = OrdinaryDiffEq.ContinuousCallback(condition, affect!, nothing) # first affect is to stop when going from neg to pos, second affect is to stop when going from pos to neg

    sol = solve(prob, Vern9(), abstol=1e-12, reltol=1e-12, callback=cb) # solve the problem
    x = sol.u[end][1:6] # final state
    t = sol.t[end] # final time
    Φ = reshape(sol.u[end][7:42],6,6) # final STM
    ẋ = CR3BPdynamicsBar(x, μ, 0)
    Dt = -(Dσ*ẋ)\Dσ*Φ
    DP = Φ + ẋ*Dt
    return x, t, DP
end

"""
    map_CR3BP_reduced(z, μ, H, tmax=10)

Poincare map on the reduced 4D space (q₁, q₂, p₁, p₂).
Returns (fz, Dfz) where fz is the mapped point and Dfz is the 4×4 Jacobian.
"""
function map_CR3BP_reduced(z, μ, H, tmax=10)
    # Lift z ∈ ℝ⁴ to x ∈ ℝ⁶
    x = γ(z, H, μ)

    # Apply Poincare map in full space
    fx, t, DPx = P(x, μ, tmax)

    # Project back to reduced space
    fz = γ⁻¹(fx)

    # Chain rule: Dfz = Dγ⁻¹(fx) * DPx * Dγ(z)
    Dfz = Dγ⁻¹(fx) * DPx * Dγ(z, μ)

    return fz, Dfz
end

"""
    shift_fourier(paramF::Array{ComplexF64, 4}, ω::Vector{Float64})

Shift Fourier coefficients by rotation vector ω.
For a function K(θ) with Fourier series K̂(k), computes K(θ+ω) by multiplying K̂(k) by exp(2πi k·ω).

# Arguments
- `paramF`: Fourier coefficients [DMAP, 1, DTOR, nelem] - complex array
- `ω`: Rotation frequencies [DTOR] - real vector

# Returns
- `KshiftF`: Shifted Fourier coefficients with same shape
"""
function shift_fourier(paramF::Array{ComplexF64, 4}, ω::Vector{Float64})
    DMAP, _, n1, n2 = size(paramF)
    nn = [n1, n2]
    DTOR = length(ω)

    KshiftF = copy(paramF)

    # Loop over all Fourier modes
    for i1 in 1:n1, i2 in 1:n2
        # Convert to signed frequencies: k ∈ [1, nn[i]] → k ∈ [-nn[i]/2, nn[i]/2]
        # Following FFT convention: frequencies 0...N/2-1, -N/2...-1
        # Note: Julia uses 1-based indexing, so we subtract 1 first
        k1 = i1 - 1  # Now in [0, n1-1]
        k2 = i2 - 1  # Now in [0, n2-1]

        # Shift to signed frequencies
        if k1 > n1 ÷ 2
            k1 -= n1
        end
        if k2 > n2 ÷ 2
            k2 -= n2
        end

        k_vec = [k1, k2]

        # Compute phase shift: exp(2πi k·ω)
        phase = exp(2π * im * dot(k_vec, ω))

        # Apply shift to all components
        for j in 1:DMAP
            KshiftF[j, 1, i1, i2] *= phase
        end
    end

    return KshiftF
end

"""
    norm_fourier(paramF::Array{ComplexF64, 4})

Compute L1 norm of Fourier coefficients (sum of absolute values).
"""
function norm_fourier(paramF::Array{ComplexF64, 4})
    return sum(abs, paramF)
end

"""
    kam_torus!(paramR, paramF, ω, μ, H; toltail=1e-12, tolinva=1e-12, Case=2)

Perform one KAM correction iteration on an invariant torus parameterization.

# Arguments
- `paramR`: Real-space torus parameterization [DMAP, 1, nn[1], nn[2]] - modified in place
- `paramF`: Fourier-space torus parameterization [DMAP, 1, nn[1], nn[2]] - modified in place
- `ω`: Rotation frequencies [ω₁, ω₂]
- `μ`: CR3BP mass parameter
- `H`: Hamiltonian value

# Keyword Arguments
- `toltail`: Tolerance for Fourier tail convergence (default: 1e-12)
- `tolinva`: Tolerance for invariance error (default: 1e-12)
- `Case`: Frame construction method (1: fixed normal, 2: metric-based) (default: 2)

# Returns
- `error`: L1 norm of invariance error
- `converged`: true if error < tolinva, false otherwise
- `tail_too_large`: true if Fourier tail exceeds toltail

# Algorithm Steps
1. STEP 0: Check Fourier tail size
2. STEP 1: Compute invariance error E = Map(K(θ)) - K(θ+ω)
3. STEP 2: Build symplectic frame (tangent L, normal N)
4. STEP 3: Solve cohomological equations for corrections
5. STEP 4: Update parameterization K_new = K + corrections

This implements Algorithm 4.32 from Haro et al. (2016).
"""
function kam_torus!(paramR, paramF, ω, μ, H;
                    toltail=1e-12, tolinva=1e-12, Case=2)

    DMAP, _, nn1, nn2 = size(paramR)
    DTOR = length(ω)
    nn = [nn1, nn2]
    nelem = prod(nn)

    @assert DMAP == 4 "Expected 4D reduced phase space"
    @assert DTOR == 2 "Expected 2D torus"

    println("#     - Size of the grid: $(nn[1]) $(nn[2])")

    # ========================================================================
    # STEP 0: Evaluation of the tail
    # ========================================================================

    # Compute Fourier tail (equation 4.104 in Haro et al. 2016)
    # Tail measures decay of high-frequency Fourier coefficients
    tg = zeros(DTOR)
    tail_too_large = false

    # For each angular direction, compute tail as sum of coefficients in middle band
    # Following C++ implementation: sum coefficients where nn[i]/4 ≤ index[i] ≤ 3*nn[i]/4
    # This measures decay of mid-frequency coefficients (not just high-frequency)
    for dim in 1:DTOR
        tail_sum = 0.0
        # Collect coefficients in the middle frequency band for this dimension
        for i1 in 1:nn[1], i2 in 1:nn[2]
            # Convert to 0-based index for comparison with C++
            idx = (dim == 1) ? (i1 - 1) : (i2 - 1)
            n = nn[dim]

            # Check if in middle band: n/4 ≤ idx ≤ 3n/4
            if n ÷ 4 <= idx && idx <= 3 * n ÷ 4
                # Sum absolute values of all components at this grid point
                for j in 1:DMAP
                    tail_sum += abs(paramF[j, 1, i1, i2])
                end
            end
        end

        # Average by number of elements (factor of 2 from C++ symmetry consideration)
        tg[dim] = 2.0 * tail_sum / nelem

        if tg[dim] > toltail
            tail_too_large = true
        end
    end

    println("#     - Tails of the parameterization: $(tg[1]) $(tg[2])")

    if tail_too_large
        println("#     - Tail too large, need finer grid resolution")
        return NaN, false, true
    end

    # Clean high-frequency modes (optional truncation)
    # In C++ this is done with clean(paramF) - we skip for now

    # ========================================================================
    # STEP 1: Evaluation of the invariance error
    # ========================================================================

    println("#     - Computing invariance error...")

    # Allocate arrays for map evaluation at each grid point
    FparamR = zeros(ComplexF64, DMAP, 1, nn[1], nn[2])  # Map(K(θ))

    # Evaluate map at each grid point
    for i1 in 1:nn[1], i2 in 1:nn[2]
        # Extract state vector z at this grid point
        z = real([paramR[j, 1, i1, i2] for j in 1:DMAP])

        # Evaluate Poincare map
        fz, Dfz = map_CR3BP_reduced(z, μ, H)

        # Store result
        for j in 1:DMAP
            FparamR[j, 1, i1, i2] = fz[j]
        end
    end

    # Compute K(θ+ω) by shifting in Fourier space
    # First convert to Fourier space
    FparamF = fft(FparamR, (3, 4))

    # Shift by ω
    KshiftF = shift_fourier(paramF, ω)

    # Convert shifted torus back to real space
    KshiftR = ifft(KshiftF, (3, 4))

    # Add angle offset to first DTOR components
    # K_shift(θ) = K(θ+ω) + (θ+ω) in the angular coordinates
    for i1 in 1:nn[1], i2 in 1:nn[2]
        # Grid angles θ[i] = (i-1)/nn[i] for i = 1, 2, ..., nn[i]
        θ1 = (i1 - 1) / nn[1]
        θ2 = (i2 - 1) / nn[2]

        # Add angle offset: θ + ω (only for first DTOR components)
        # Note: C++ line 634-635 was commented out, but line 662 adds it back
        for j in 1:DTOR
            θ_j = j == 1 ? θ1 : θ2
            KshiftR[j, 1, i1, i2] += θ_j + ω[j]
        end
    end

    # Compute invariance error: E = Map(K(θ)) - K(θ+ω)
    ErrorR = FparamR - KshiftR

    # Convert to Fourier space and compute norm
    ErrorF = fft(ErrorR, (3, 4))
    error = norm_fourier(ErrorF)

    println("#     - Error of invariance: $error")

    if error < tolinva
        println("#     - No correction is needed!")
        return error, true, false
    end

    # ========================================================================
    # STEP 2: Build symplectic frame (TODO)
    # ========================================================================

    println("#     - Building symplectic frame... (TODO)")

    # TODO: Implement tangent frame L = ∂K/∂θ
    # TODO: Implement normal frame N

    # ========================================================================
    # STEP 3: Solve cohomological equations (TODO)
    # ========================================================================

    println("#     - Solving cohomological equations... (TODO)")

    # TODO: Implement small divisor solver
    # TODO: Implement twist matrix inversion for average component

    # ========================================================================
    # STEP 4: Update parameterization (TODO)
    # ========================================================================

    println("#     - Updating parameterization... (TODO)")

    # TODO: Apply corrections K_new = K + L·ξ_L + N·ξ_N

    return error, false, false
end
