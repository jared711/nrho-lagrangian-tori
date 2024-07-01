using LinearAlgebra

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
    @assert isapprox(x[3], 0) "x[3] must be zero."
    # @assert isapprox(computeH(x, μ) ≈ H) "The Hamiltonian is not conserved."
    z = x[[1,2,4,5]]
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
    
    sol = solve(prob, TsitPap8(), abstol=1e-12, reltol=1e-12, callback=cb) # solve the problem
    x = sol.u[end][1:6] # final state
    t = sol.t[end] # final time
    Φ = reshape(sol.u[end][7:42],6,6) # final STM
    ẋ = CR3BPdynamicsBar(x, μ, 0)
    Dt = -(Dσ*ẋ)\Dσ*Φ
    DP = Φ + ẋ*Dt
    return x, t, DP
    # sys_Barcelona = System(ENCELADUS, SATURN, ENCELADUS.a, ENCELADUS.T,
    # ENCELADUS.m*G, SATURN.m*G, SATURN.m/(ENCELADUS.m+SATURN.m), ENCELADUS.R, ENCELADUS.R,
    # ENCELADUS.m*G, SATURN.m*G, SATURN.m/(ENCELADUS.m+SATURN.m), ENCELADUS.R, ENCELADUS.R,
    # ENCELADUS.m+SATURN.m, ENCELADUS.a, ENCELADUS.T/2π, ENCELADUS.a/(ENCELADUS.T/2π), ENCELADUS.a/(ENCELADUS.T/2π)^2, string(ENCELADUS.name,"/",SATURN.name))

    # plot(sol,vars=(1,2),label="Halo Orbit",legend=false,xaxis="x",yaxis="y",zaxis="z",color=:blue) # plot the halo orbit in the x-y plane
    # plot!(sys_Barcelona,sec=false,Lpts=false, lims=:auto)
end
