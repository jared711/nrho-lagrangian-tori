"""
    CR3BPdynamicsBar(x, μ, t)

Compute time derivative of state vector `x = [q; p]` {NON, NON} in the rotating frame of
the normalized CR3BP where `μ` is the CR3BP mass parameter μ₂/(μ₁+μ₂) {NON} and `t` is time
{NON}. This uses the Barcelona equations. Note that the secondary body is on the left
"""
function CR3BPdynamicsBar(x, μ, t) #Three body dynamics function in Hamiltonian form
    q₁, q₂, q₃, p₁, p₂, p₃ = x
    r₁³= ((q₁ - μ)^2     + q₂^2 + q₃^2)^1.5 # distance to m1, LARGER MASS
    r₂³= ((q₁ + 1 - μ)^2 + q₂^2 + q₃^2)^1.5 # distance to m2, smaller mass

    xdot = zeros(6)
    xdot[1] = p₁ + q₂
    xdot[2] = p₂ - q₁
    xdot[3] = p₃
    xdot[4] = -((1 - μ)*(q₁ - μ)/r₁³) - (μ*(q₁ + 1 - μ)/r₂³) + p₂
    xdot[5] = -((1 - μ)*q₂      /r₁³) - (μ*q₂          /r₂³) - p₁
    xdot[6] = -((1 - μ)*q₃      /r₁³) - (μ*q₃          /r₂³)
    
    return xdot
end

# """
#     CR3BPdynamicsBar(x, sys::System, t)

# Compute time derivative of state vector `x = [q; p]` {NON, NON} in the rotating frame of
# the normalized CR3BP where `sys` is the CR3BP system and `t` is time {NON}. Uses the Barcelona convention
# """
# function CR3BPdynamicsBar(x, sys::System, t) #Three body dynamics function in Hamiltonian form
#     return CR3BPdynamicsBar(x, sys.μ, t)
# end

"""
    CR3BPdynamicsBar!(rvdot, rv, μ, t)

In-place version of `CR3BPdynamicsBar(rv, μ, t)`.
"""
function CR3BPdynamicsBar!(rvdot,rv,μ,t) #Three body dynamics in Earth/Moon System
    rvdot[:] = CR3BPdynamicsBar(rv,μ,t)
    return nothing
end

# """
#     CR3BPdynamicsBar!(rvdot, rv, sys::System, t)

# In-place version of `CR3BPdynamicsBar(rv, sys::System, t)`.
# """
# function CR3BPdynamicsBar!(rvdot,rv,sys::System,t) #Three body dynamics in Earth/Moon System
#     rvdot[:] = CR3BPdynamicsBar(rv,sys,t)
#     return nothing
# end

"""
    CR3BPjacBar(x, μ)

Compute Jacobian of time derivative w.r.t. state vector [6x6]
"""
function CR3BPjacBar(x, μ)
    q₁, q₂, q₃ = x

    r₁³= ((q₁-μ)^2   + q₂^2 + q₃^2)^1.5 # distance to m1, LARGER MASS
    r₂³= ((q₁+1-μ)^2 + q₂^2 + q₃^2)^1.5 # distance to m2, smaller mass
    r₁⁵= ((q₁-μ)^2   + q₂^2 + q₃^2)^2.5 # distance to m1, LARGER MASS
    r₂⁵= ((q₁+1-μ)^2 + q₂^2 + q₃^2)^2.5 # distance to m2, smaller mass

    ∂ṗ₁_∂q₁ = (1-μ)*(3*(q₁-μ)^2/r₁⁵ - 1/r₁³) + μ*(3*(q₁+1-μ)^2/r₂⁵ - 1/r₂³)
    ∂ṗ₂_∂q₂ = (1-μ)*(3*    q₂^2/r₁⁵ - 1/r₁³) + μ*(3*      q₂^2/r₂⁵ - 1/r₂³)
    ∂ṗ₃_∂q₃ = (1-μ)*(3*    q₃^2/r₁⁵ - 1/r₁³) + μ*(3*      q₃^2/r₂⁵ - 1/r₂³)
    
    ∂ṗ₁_∂q₂ = 3*(1-μ)*(q₁-μ)*q₂/r₁⁵ + 3*μ*(q₁+1-μ)*q₂/r₂⁵
    ∂ṗ₁_∂q₃ = 3*(1-μ)*(q₁-μ)*q₃/r₁⁵ + 3*μ*(q₁+1-μ)*q₃/r₂⁵
    ∂ṗ₂_∂q₃ = 3*(1-μ)*    q₂*q₃/r₁⁵ + 3*μ*      q₂*q₃/r₂⁵


    F = [      0        1        0     1     0     0 ;
              -1        0        0     0     1     0 ;
               0        0        0     0     0     1 ;
         ∂ṗ₁_∂q₁  ∂ṗ₁_∂q₂  ∂ṗ₁_∂q₃     0     1     0 ;
         ∂ṗ₁_∂q₂  ∂ṗ₂_∂q₂  ∂ṗ₂_∂q₃    -1     0     0 ;
         ∂ṗ₁_∂q₃  ∂ṗ₂_∂q₃  ∂ṗ₃_∂q₃     0     0     0 ]
    return F
end

# """
#     CR3BPjacBar(rv, sys)

# Compute Jacobian of time derivative w.r.t. state vector [6x6]
# """
# CR3BPjacBar(rv, sys::System) = CR3BPjacBar(rv, sys.μ)

"""
    CR3BPstmBar(w, μ, t)

Compute time derivative of state vector `w = [q; p; vec(Φ)]` {NON; NON; NON} in the rotating
frame of the normalized CR3BP. `vec(Φ)` is the vectorized state transition matrix while `μ`
is the CR3BP mass parameter μ₂/(μ₁+μ₂) {NON} and `t` is time {NON}. Using the Barcelona convention
"""
function CR3BPstmBar(w, μ, t) #Three body dynamics in Earth/Moon System
    x = w[1:6]
    Φ = reshape(w[7:42],6,6)

    F = CR3BPjacBar(x, μ)

    Φdot = F*Φ
    wdot = zeros(42)
    wdot[1:6] = CR3BPdynamicsBar(x,μ,t)
    wdot[7:42] = reshape(Φdot, 36, 1)
    return wdot
end

# """
#     CR3BPstmBar(w, sys, t)

# Compute time derivative of state vector `w = [q; p; vec(Φ)]` {NON; NON; NON} in the rotating
# frame of the normalized CR3BP. `vec(Φ)` is the vectorized state transition matrix while
# `sys` is the CR3BP system and `t` is time {NON}.
# """
# function CR3BPstmBar(w,sys::System,t)
#     return CR3BPstmBar(w,sys.μ,t)
# end

"""
    CR3BPstmBar!(wdot, w, μ, t)

In-place version of `CR3BPstmBar(w, μ, t)`.
"""
function CR3BPstmBar!(wdot,w,μ,t)
    wdot[:] = CR3BPstmBar(w,μ,t)
    return nothing
end

# """
#     CR3BPstmBar!(wdot, w, sys::System, t)

# In-place version of `CR3BPstmBar(w, sys::System, t)`.
# """
# function CR3BPstmBar!(wdot,w,sys::System,t) #Three body dynamics in Earth/Moon System
#     wdot[:] = CR3BPstmBar(w,sys,t)
#     return nothing
# end


