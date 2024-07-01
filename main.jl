using ThreeBodyProblem
using OrdinaryDiffEq
using Plots
using CSV
using DataFrames
using LinearAlgebra

include("util.jl")
include("dynamics.jl")
include("parameterization-method.jl")

# Read in the initial conditions to a dataframe
# initial_conditions = read_initial_conditions("initial_conditions.txt")
NRHO_L2 = CSV.read("NRHO_L2.csv", DataFrame, normalizenames=true)

# Take the initial conditions of the row defined by idx
idx = 72
HaloIC = NRHO_L2[NRHO_L2.:Id.== idx, :]
# rv₀ = Vector(NRHO_L2[idx, ["x0_LU_", "y0_LU_", "z0_LU_", "vx0_LU_TU_", "vy0_LU_TU_", "vz0_LU_TU_"]])
rv₀ = Vector(HaloIC[1,["x0_LU_", "y0_LU_", "z0_LU_", "vx0_LU_TU_", "vy0_LU_TU_", "vz0_LU_TU_"]])
ν = HaloIC[1, "Stability_index"]
μ = HaloIC[1, "Mass_ratio"]
C₀ = HaloIC[1, "Jacobi_constant_LU2_TU2_"]
T₀ = HaloIC[1, "Period_TU_"]
H₀ = -C₀/2

# convert to the Hamiltonian coordinates (Barcelona Convention)
x₀ = rv2pq(rv₀) # initial condition of halo orbit (normally at apoapsis)
@assert isapprox(computeH(x₀, μ), H₀, atol=1e-12) # double check that the Hamiltonian and Jacobi Constant Agree

# Take the Poincare Map
# xₚ, tₚ = P(x₀, μ)

# We'll integrate the halo orbit and its state transition matrix (STM)
Φ₀ = I(6) # Initialization of the STM, Φ₀ = I
w₀ = [x₀; reshape(Φ₀,36,1)] # Reshape the matrix into a vector and append it to the state vector
tspan = (0.,T₀) # integrate from 0 to T₀
prob_halo = ODEProblem(CR3BPstmBar!,w₀,tspan,μ) # CR3BPstm! is our in-place dynamics function for state and STM
halo = solve(prob_halo,TsitPap8(),abstol=1e-12,reltol=1e-12) # solve the problem

# # # We'll integrate the halo orbit and its state transition matrix (STM)
# Φ₀ = I(6) # Initialization of the STM, Φ₀ = I
# w₀ = [rv₀; reshape(Φ₀,36,1)] # Reshape the matrix into a vector and append it to the state vector
# tspan = (0.,T₀) # integrate from 0 to T₀
# prob_halo = ODEProblem(CR3BPstm!,w₀,tspan,μ) # CR3BPstm! is our in-place dynamics function for state and STM
# halo = solve(prob_halo,TsitPap8(),abstol=1e-12,reltol=1e-12) # solve the problem

# write trajectory to CSV file
halo_df = DataFrame(halo, [:t, :q1, :q2, :q3, :p1, :p2, :p3, :Phi11, :Phi12, :Phi13, :Phi14, :Phi15, :Phi16, :Phi21, :Phi22, :Phi23, :Phi24, :Phi25, :Phi26, :Phi31, :Phi32, :Phi33, :Phi34, :Phi35, :Phi36, :Phi41, :Phi42, :Phi43, :Phi44, :Phi45, :Phi46, :Phi51, :Phi52, :Phi53, :Phi54, :Phi55, :Phi56, :Phi61, :Phi62, :Phi63, :Phi64, :Phi65, :Phi66])
CSV.write("halo.csv", halo_df, writeheader=true, delim=' ')

w = halo.u[end] # final state vector
x = w[1:6] # final state
Φ = reshape(w[7:42],6,6) # final STM


# import ThreeBodyProblem.eig
λ, V = eig(Φ)
@assert isapprox(stability_index(Φ), ν, atol=1e-6)

# Time to approximate an invariant 2 torus on the stroboscopic map
# We'll use the eigenvectors of the STM to approximate the torus
# N needs to be even!
N = 32 # Number of points on the invariant circle (THIS SHOULD BE A POWER OF 2)
# N is also the number of frequencies that we will break our u function into
n = 6 # Number of dimensions of the system

eig_idx1 = 4
eig_idx2 = 6

ρ₁ = atan(imag(λ[eig_idx1]),real(λ[eig_idx1])) # Angle of the first eigenvector
ρ₂ = atan(imag(λ[eig_idx2]),real(λ[eig_idx2])) # Angle of the second eigenvector
ω₁ = ρ₁/T₀
ω₂ = ρ₂/T₀

θ = 2π*(0:N-1)/N # Angles for the invariant circle
α = 1e-5 # parameter to control the size of the invariant circle
β = 1e-5
u = zeros(N^2,n) # Initialize the invariant circle
uP = zeros(N^2,n) # Initialize the invariant circle on the Poincare Section
idx1_vec = []
idx2_vec = []
uidx = 0
@time for (idx1, θ₁) in enumerate(θ)
    for (idx2, θ₂) in enumerate(θ)
        uidx += 1
        u[uidx,:] = α*(cos(θ₁)*real(V[:,eig_idx1]) - sin(θ₁)*imag(V[:,eig_idx1])) + β*(cos(θ₂)*real(V[:,eig_idx2]) - sin(θ₂)*imag(V[:,eig_idx2])) # Initial guess for the invariant circle
        u[uidx,:] += x
        push!(idx1_vec, idx1)
        push!(idx2_vec, idx2)

        uP[uidx,:],_,_ = P(u[uidx,:],μ)
    end
end

H = [computeH(uP[i,:],μ) for i in 1:N^2] # Compute the Jacobi constant of the invariant circle
paramsdf = DataFrame(
    toltail = 1.0e-11, # was 1.e-13
    tolinva = 1.0e-12,
    tolinte = 1.0e-17,
    ω₁ = ω₁,
    ω₂ = ω₂,
    ε = 0.0,
    λ₁ = H₀,
    λ₂ = μ,
    n₁ = N,
    n₂ = N,
    dε = 5e-5,
    auxint = 1
)

uPdf = DataFrame(uP[:,[1,2,4,5]],["q₁","q₂","p₁","p₂"]) # Display the invariant circle in the Hamiltonian coordinates on the Poincare Map
uPdf.idx1 = idx1_vec
uPdf.idx2 = idx2_vec
select!(uPdf, "idx1", "idx2", :)
CSV.write("approxQPO.csv", paramsdf, delim=' ',writeheader=false, append=false)
CSV.write("approxQPO.csv", uPdf, writeheader=false, delim=' ',append=true)


udf = DataFrame(u[:,[1,2,4,5]],["q₁","q₂","p₁","p₂"]) # Display the invariant circle in the Hamiltonian coordinates
udf.idx1 = idx1_vec
udf.idx2 = idx2_vec
select!(udf, "idx1", "idx2", :)
CSV.write("approxQPO.csv", paramsdf, delim=' ',writeheader=false, append=false)
CSV.write("approxQPO.csv", udf, writeheader=false, delim=' ',append=true)

### Stuff from the kam_torus function



### Plotting stuff to verify that the Pmap is working correctly
pmap_plot = scatter(udf.q₁, udf.q₂, xlabel="q₁", ylabel="q₂", label="2D torus approx", legend=true, title="Invariant Circle in Hamiltonian Coordinates", markercolor=:blue, marker=:x)
plot!(pmap_plot,[x[1]],[x[2]],label="x₀",shape=:o,markercolor=:red,legend=true,line=false) # Plot the initial condition

sys_Barcelona = System(ENCELADUS, SATURN, ENCELADUS.a, ENCELADUS.T,
ENCELADUS.m*ThreeBodyProblem.G, SATURN.m*ThreeBodyProblem.G, SATURN.m/(ENCELADUS.m+SATURN.m), ENCELADUS.R, ENCELADUS.R,
ENCELADUS.m*ThreeBodyProblem.G, SATURN.m*ThreeBodyProblem.G, SATURN.m/(ENCELADUS.m+SATURN.m), ENCELADUS.R, ENCELADUS.R,
ENCELADUS.m+SATURN.m, ENCELADUS.a, ENCELADUS.T/2π, ENCELADUS.a/(ENCELADUS.T/2π), ENCELADUS.a/(ENCELADUS.T/2π)^2, string(ENCELADUS.name,"/",SATURN.name))

plot_haloWithTorus = plot(halo,vars=(1,2),label="Halo Orbit",legend=false,xaxis="x",yaxis="y",zaxis="z",color=:blue) # plot the halo orbit in the x-y plane
plot!(plot_haloWithTorus,sys_Barcelona,sec=false,Lpts=false, lims=:auto)
scatter!(plot_haloWithTorus,udf.q₁, udf.q₂, xlabel="q₁", ylabel="q₂", label="2D torus approx", legend=true, title="Invariant Circle in Hamiltonian Coordinates", markercolor=:red, marker=:x)
plot!(plot_haloWithTorus,[x[1]],[x[2]],label="x₀",shape=:o,markercolor=:blue,legend=true,line=false) # Plot the initial condition

u_pmap = [P(u[i,:],μ) for i = 1:N^2]
u_pmap_x = [u_pmap[i][1] for i = 1:N^2]
u_pmap_y = [u_pmap[i][2] for i = 1:N^2]

scatter( u_pmap_x, u_pmap_y,label="2D torus approx",legend=true,markercolor=:red,marker=:x) # Plot the invariant circle



T₀ = [γ⁻¹(x₀ + u[:,i]) for i in 1:N^2] # Initial guess for the invariant circle in the Hamiltonian coordinates
CSV.write("T₀.csv", DataFrame(T₀), writeheader=false)

plot_u = plot(u, xlabel="X [NON]",ylabel="Y [NON]", zlabel= "Z [NON]", legend=true,label="u",title="Approximate Invariant Circle",linecolor=:blue, marker=:x); # Plot the invariant circle
scatter!(plot_u, [u[1][1]],[u[1][2]],[u[1][3]],label="u[1]",shape=:o,markercolor=:blue) # Plot an "x" on the first point of the invariant circle
     


z₀ = γ⁻¹(x₀)


# Plot the halo orbit
pxy = plot(halo,idxs=(1,2),label="Halo Orbit",legend=false,xaxis="x",yaxis="y"); # plot the halo orbit in the x-y plane
plot!(System(ENCELADUS,SATURN),sec=false,Lpts=false, lims=:auto);

pyz = plot(halo,idxs=(2,3),label="Halo Orbit",legend=false,xaxis="y",yaxis="z"); # plot the halo orbit in the y-z plane
plot!(System(ENCELADUS,SATURN),sec=false,Lpts=false, lims=:auto, center=[0,0,0]);

pxz = plot(halo,idxs=(1,3),label="Halo Orbit", legend=false,xaxis="x",yaxis="z"); # plot the halo orbit in the x-z plane
plot!(System(ENCELADUS,SATURN),sec=false,Lpts=false, lims=:auto);

# pall = plot(halo,idxs=(1,2,3),legend=false,title="Halo Orbit",label="Halo Orbit"); # plot the halo orbit in 3D
# plot!(pall,sys,planar=false,prim=false,Lpts=false,lims=:auto);

plot_halo = plot(pxy,pyz,pxz,layout=(1,3),legend=:outertop) # plot all of the plots in a 2x2 grid
# plot_halo = plot(pxy,pyz,pxz,pall,layout=(1,4),title="Halo Orbit") # plot all of the plots in a 2x2 grid
