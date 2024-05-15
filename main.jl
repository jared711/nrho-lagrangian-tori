using ThreeBodyProblem
using DifferentialEquations
using Plots
using CSV
using DataFrames
using LinearAlgebra

sys = saturn_enceladus() # Define the system we're working with

# Read in the initial conditions to a dataframe
initial_conditions = read_initial_conditions("initial_conditions.txt")
NRHO_L1 = CSV.read("NRHO_L1.csv", DataFrame)

# Take the initial conditions of the row defined by idx
idx = 1
rv₀ = Vector(NRHO_L1[idx, ["x0 (LU) ", "y0 (LU) ", "z0 (LU) ", "vx0 (LU/TU) ", "vy0 (LU/TU) ", "vz0 (LU/TU) "]])
C₀ = NRHO_L1[idx, "Jacobi constant (LU2/TU2) "]
T₀ = NRHO_L1[idx, "Period (TU) "]

# We'll integrate the halo orbit and its state transition matrix (STM)
Φ₀ = I(6) # Initialization of the STM, Φ₀ = I
w₀ = [rv₀; reshape(Φ₀,36,1)] # Reshape the matrix into a vector and append it to the state vector
tspan = (0.,T₀) # integrate from 0 to T₀
prob_halo = ODEProblem(CR3BPstm!,w₀,tspan,sys) # CR3BPstm! is our in-place dynamics function for state and STM
halo = solve(prob_halo,TsitPap8(),abstol=1e-12,reltol=1e-12) # solve the problem

# Plot the halo orbit
pxy = plot(halo,idxs=(1,2),label="Halo Orbit",legend=false,xaxis="x",yaxis="y"); # plot the halo orbit in the x-y plane
plot!(sys,prim=false,Lpts=false, lims=:auto);

pyz = plot(halo,idxs=(2,3),label="Halo Orbit",legend=false,xaxis="y",yaxis="z"); # plot the halo orbit in the y-z plane
plot!(sys,prim=false,Lpts=false, lims=:auto, center=[0,0,0]);

pxz = plot(halo,idxs=(1,3),label="Halo Orbit", legend=false,xaxis="x",yaxis="z"); # plot the halo orbit in the x-z plane
plot!(sys,prim=false,Lpts=false, lims=:auto);

# pall = plot(halo,idxs=(1,2,3),legend=false,title="Halo Orbit",label="Halo Orbit"); # plot the halo orbit in 3D
# plot!(pall,sys,planar=false,prim=false,Lpts=false,lims=:auto);

plot_halo = plot(pxy,pyz,pxz,layout=(1,3),legend=:outertop) # plot all of the plots in a 2x2 grid
# plot_halo = plot(pxy,pyz,pxz,pall,layout=(1,4),title="Halo Orbit") # plot all of the plots in a 2x2 grid

