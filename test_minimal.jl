using Pkg
Pkg.activate(".")

using CSV, DataFrames, FFTW, LinearAlgebra, OrdinaryDiffEq, Statistics

include("src/julia/util.jl")
include("src/julia/dynamics.jl")
include("src/julia/parameterization-method.jl")

println("===== Minimal KAM Torus Test =====\n")

# Load test data
datafile = "data/initial_conditions/approxQPO.csv"
params_line = readline(datafile)
params = parse.(Float64, split(params_line))

μ = params[8]
H = params[7]
ω = [params[4], params[5]]
n1, n2 = Int(params[9]), Int(params[10])

println("Parameters:")
println("  μ = $μ")
println("  H = $H")
println("  ω = $ω")
println("  Grid: $n1 × $n2\n")

# Load grid data
df = CSV.read(datafile, DataFrame, header=false, skipto=2, delim=' ')
println("Loaded $(nrow(df)) grid points\n")

# Test single Poincare map evaluation
println("Testing single Poincare map evaluation...")
z_test = [df[1, 3], df[1, 4], df[1, 5], df[1, 6]]  # First grid point
println("  Input z = $z_test")

try
    fz, Dfz = map_CR3BP_reduced(z_test, μ, H, 10.0)
    println("  Output fz = $fz")
    println("  Jacobian size = $(size(Dfz))")
    println("  ✓ Poincare map works!\n")
catch e
    println("  ✗ ERROR: $e\n")
    rethrow(e)
end

# Test Fourier operations
println("Testing Fourier operations...")
test_array = zeros(ComplexF64, 4, 1, n1, n2)
test_array[1, 1, 1, 1] = 1.0 + 0.0im  # DC component
println("  Test array created")

test_fft = fft(test_array, (3, 4))
println("  FFT computed")

test_shift = shift_fourier(test_fft, ω)
println("  Shift computed")

test_norm = norm_fourier(test_shift)
println("  Norm = $test_norm")
println("  ✓ Fourier operations work!\n")

println("===== All basic operations successful! =====")
