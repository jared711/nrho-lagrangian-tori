using CSV
using DataFrames
using FFTW
using LinearAlgebra
using OrdinaryDiffEq
using Statistics

include("util.jl")
include("dynamics.jl")
include("parameterization-method.jl")

"""
Load test data from approxQPO.csv and run one KAM torus iteration.
"""
function test_kam_torus()
    println("=" ^ 70)
    println("Testing KAM Torus Algorithm - STEP 1 Implementation")
    println("=" ^ 70)

    # Load test data
    datafile = "data/initial_conditions/approxQPO.csv"
    println("\nLoading test data from: $datafile")

    # Read parameters from first line
    params_line = readline(datafile)
    params = parse.(Float64, split(params_line))

    toltail = params[1]
    tolinva = params[2]
    tolinte = params[3]
    ω1 = params[4]
    ω2 = params[5]
    ε = params[6]
    H = params[7]  # This is λ₁ (Hamiltonian)
    μ = params[8]  # This is λ₂ (mass parameter)
    n1 = Int(params[9])
    n2 = Int(params[10])
    dε = params[11]
    auxint = Int(params[12])

    ω = [ω1, ω2]

    println("\nParameters:")
    println("  toltail = $toltail")
    println("  tolinva = $tolinva")
    println("  tolinte = $tolinte")
    println("  ω = [$ω1, $ω2]")
    println("  H = $H")
    println("  μ = $μ")
    println("  Grid dimensions: $n1 × $n2")

    # Load torus data (skip first line which is parameters)
    df = CSV.read(datafile, DataFrame, header=false, skipto=2, delim=' ')

    # Data format: idx1 idx2 q₁ q₂ p₁ p₂
    println("\nLoaded $(nrow(df)) grid points")
    println("Sample data (first 3 rows):")
    println(first(df, 3))

    # Verify grid dimensions
    @assert nrow(df) == n1 * n2 "Expected $(n1*n2) points, got $(nrow(df))"

    # Convert to parameterization arrays
    # paramR[component, 1, i1, i2] where component ∈ {q₁, q₂, p₁, p₂}
    DMAP = 4  # 4D reduced phase space
    DTOR = 2  # 2D torus

    paramR = zeros(ComplexF64, DMAP, 1, n1, n2)

    for row in eachrow(df)
        i1 = Int(row[1])  # idx1 (1-indexed)
        i2 = Int(row[2])  # idx2 (1-indexed)
        q1 = row[3]
        q2 = row[4]
        p1 = row[5]
        p2 = row[6]

        paramR[1, 1, i1, i2] = q1
        paramR[2, 1, i1, i2] = q2
        paramR[3, 1, i1, i2] = p1
        paramR[4, 1, i1, i2] = p2
    end

    println("\nConverting to Fourier space...")
    # Convert to Fourier representation
    # FFT over dimensions 3 and 4 (the angular grid dimensions)
    paramF = fft(paramR, (3, 4))

    println("paramR statistics:")
    println("  min = $(minimum(real.(paramR)))")
    println("  max = $(maximum(real.(paramR)))")
    println("  mean = $(mean(real.(paramR)))")

    println("\nparamF statistics:")
    println("  |paramF|₀₀ (DC component) = $(abs(paramF[1,1,1,1]))")
    println("  max |paramF| = $(maximum(abs.(paramF)))")

    # Run KAM torus iteration (STEP 1 only for now)
    println("\n" * "=" ^ 70)
    println("Running KAM Torus Iteration (STEP 1 - Invariance Error)")
    println("=" ^ 70)

    error, converged, tail_too_large = kam_torus!(
        paramR, paramF, ω, μ, H;
        toltail=toltail,
        tolinva=tolinva,
        Case=2
    )

    println("\n" * "=" ^ 70)
    println("Results:")
    println("=" ^ 70)
    println("  Invariance error: $error")
    println("  Converged: $converged")
    println("  Tail too large: $tail_too_large")

    # Expected from C++ output: error ≈ 2.593
    println("\n  Expected error (from C++): ≈ 2.593")
    if !isnan(error)
        rel_error = abs(error - 2.593) / 2.593
        println("  Relative difference: $(100*rel_error)%")

        if rel_error < 0.1
            println("\n✓ PASS: Error matches C++ implementation within 10%")
        else
            println("\n✗ FAIL: Error differs significantly from C++ implementation")
        end
    end

    return error, converged, tail_too_large, paramR, paramF
end

# Run the test
if abspath(PROGRAM_FILE) == @__FILE__
    error, converged, tail_too_large, paramR, paramF = test_kam_torus()
end
