# Julia Implementation of KAM Torus Algorithm

This directory contains a Julia port of the parameterization method for computing invariant tori around Near-Rectilinear Halo Orbits (NRHOs) in the Circular Restricted Three-Body Problem (CR3BP).

## Overview

The parameterization method (Haro & Luque, 2016) computes an invariant 2-torus K(θ) in 4D reduced phase space by iteratively correcting an initial approximation until the invariance equation is satisfied:

```
K(θ + ω) = Map(K(θ))
```

where:
- θ = (θ₁, θ₂) are angle variables on the 2-torus
- ω = (ω₁, ω₂) are rotation frequencies  
- Map is the Poincare map on the section q₃ = 0

## Files

### Core Implementation
- **`parameterization-method.jl`** - Main KAM algorithm
  - `kam_torus!()` - One KAM correction iteration (Algorithm 4.32)
  - `P()` - Poincare map with Jacobian computation
  - `map_CR3BP_reduced()` - 4D reduced Poincare map
  - Fourier operations: `shift_fourier()`, `norm_fourier()`
  - Coordinate transforms: `γ()`, `γ⁻¹()` between 4D and 6D

- **`dynamics.jl`** - CR3BP equations of motion
  - `CR3BPdynamicsBar()` - Hamiltonian dynamics (Barcelona convention)
  - `CR3BPjacBar()` - Jacobian matrix
  - `CR3BPstmBar!()` - Dynamics with State Transition Matrix (STM)

- **`util.jl`** - Utility functions
  - `rv2pq()` / `pq2rv()` - Convert between inertial and rotating frames
  - `computeH()` - Hamiltonian evaluation

### Testing
- **`test_kam_torus.jl`** - Test script for KAM algorithm
  - Loads test data from `data/initial_conditions/approxQPO.csv`
  - Runs one iteration and compares error with C++ implementation
  - Validates Fourier operations and map evaluation

- **`main.jl`** - Initial condition generator
  - Computes NRHO from database
  - Creates approximate torus from eigenvectors
  - Generates test data files

### Running the Tests
- **`run_test.jl`** - Convenience script (in repo root)
  ```bash
  julia --project=. run_test.jl
  ```

## Algorithm Structure

### `kam_torus!(paramR, paramF, ω, μ, H; ...)`

One KAM correction iteration following Algorithm 4.32 from Haro et al. (2016):

#### **STEP 0: Tail Evaluation**
- Computes Fourier tail size `tg[i]` for each angular direction
- If `tg[i] > toltail`, grid resolution insufficient → need more modes
- Cleans high-frequency modes to prevent numerical noise

#### **STEP 1: Invariance Error** ✅ IMPLEMENTED
- Evaluates Map(K(θ)) at all 1024 grid points
- Computes K(θ+ω) via Fourier shift
- Error: E(θ) = Map(K(θ)) - K(θ+ω)
- Returns L1 norm: `error = ||E||₁`

#### **STEP 2: Symplectic Frame** 🚧 TODO
- **Tangent frame**: L = ∂K/∂θ (Fourier differentiation)
- **Normal frame** N (two methods):
  - Case 1: User-provided transversal field
  - Case 2: Metric-based (symplectic complement via G)

#### **STEP 3: Cohomological Equations** 🚧 TODO
- Project error onto tangent/normal spaces
- Solve small divisor equations: ξ(k) = η(k) / (2πi k·ω)
- Handle average component (k=0) via twist matrix inversion

#### **STEP 4: Update Parameterization** 🚧 TODO
- Apply corrections: K_new = K + L·ξ_L + N·ξ_N
- Update paramR and paramF arrays

## Data Structures

### Torus Parameterization
```julia
paramR::Array{ComplexF64, 4}  # Real-space: [DMAP, 1, nn[1], nn[2]]
paramF::Array{ComplexF64, 4}  # Fourier-space: [DMAP, 1, nn[1], nn[2]]
```

Dimensions:
- `DMAP = 4`: Phase space components (q₁, q₂, p₁, p₂)
- `nn = [n1, n2]`: Grid resolution (typically 32×32 or 64×64)

Grid points:
- θᵢⱼ = [(i-1)/n1, (j-1)/n2] for i=1...n1, j=1...n2
- Point (i,j) stores state K(θᵢⱼ)

### Fourier Conventions
- **Forward FFT**: `paramF = fft(paramR, (3,4))`
- **Inverse FFT**: `paramR = ifft(paramF, (3,4))`
- **Frequencies**: k ∈ [0, n/2-1, -n/2, ..., -1] (standard FFTW)
- **Shift**: K(θ+ω) ↔ K̂(k) · exp(2πi k·ω)

## Dependencies

From `Project.toml`:
```toml
CSV                  # Read initial conditions
DataFrames           # Tabular data
FFTW                 # Fast Fourier Transform
LinearAlgebra        # Matrix operations
OrdinaryDiffEq       # High-precision ODE integration
ThreeBodyProblem     # CR3BP utilities
```

## Usage Example

```julia
using FFTW, LinearAlgebra, OrdinaryDiffEq

include("src/julia/util.jl")
include("src/julia/dynamics.jl")
include("src/julia/parameterization-method.jl")

# Load test data
paramR, paramF, ω, μ, H = load_test_data("data/initial_conditions/approxQPO.csv")

# Run one KAM iteration
error, converged, tail_too_large = kam_torus!(
    paramR, paramF, ω, μ, H;
    toltail = 1e-12,
    tolinva = 1e-12,
    Case = 2
)

println("Invariance error: $error")
println("Converged: $converged")
```

## Performance Notes

### Bottlenecks (from profiling)
1. **Map evaluation**: 90%+ of compute time
   - 1024 Poincare maps per iteration
   - Each requires ODE integration (50-200 steps)
   - **Optimization**: Parallelize with `Threads.@threads`

2. **FFT operations**: <5% of compute time
   - Fast even for 128×128 grids
   - Uses highly optimized FFTW library

3. **Linear algebra**: <5% of compute time
   - Mostly in cohomological solver (STEP 3)

### Optimization Strategies
- **Parallel map evaluation**: Use all CPU cores
  ```julia
  Threads.@threads for i1 in 1:nn[1]
      for i2 in 1:nn[2]
          # Evaluate map at (i1, i2)
      end
  end
  ```

- **FFTW tuning**: Create optimized plans
  ```julia
  plan = plan_fft(paramR, (3,4), flags=FFTW.MEASURE)
  paramF = plan * paramR
  ```

- **GPU acceleration**: Use CUDA.jl for large grids (128×128+)

## Test Case

**File**: `data/initial_conditions/approxQPO.csv`

**System**: Saturn-Enceladus (μ = 1.901×10⁻⁷)

**Parameters**:
```
Grid size:    32 × 32 (1024 points)
Frequencies:  ω = [0.9655, 0.4884]
Hamiltonian:  H = -1.5000
Tolerances:   toltail = 1e-11, tolinva = 1e-12
```

**Expected Results** (from C++ implementation):
- Invariance error (iteration 1): ~2.593
- Convergence after ~5-10 iterations
- Final error: <1e-12

## References

1. **Haro, A., Canadell, M., Figueras, J-Ll., Luque, A., & Mondelo, J-M.** (2016).  
   *The Parameterization Method for Invariant Manifolds: From Rigorous Results to Effective Computations*.  
   Applied Mathematical Sciences 195, Springer.  
   - Algorithm 4.32 (page ~150): KAM torus correction
   - Equation 4.104 (page ~145): Fourier tail definition

2. **C++ Reference Implementation**:  
   `src/cpp/param.cc` lines 492-834  
   See `docs/analysis/KAM_TORUS_BUG_ANALYSIS.md` for line-by-line explanation

3. **Julia Packages**:
   - [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) - High-order ODE solvers
   - [FFTW.jl](https://github.com/JuliaMath/FFTW.jl) - Fast Fourier Transform wrapper
   - [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) - Umbrella package

## Known Issues

### Current Limitations
- Only STEP 1 (invariance error) implemented
- STEP 2-4 TODO (symplectic frame, cohomological eqs, update)
- No parallelization yet (single-threaded)
- Limited error handling

### Differences from C++
1. **Array indexing**: Julia is 1-based (i=1...n), C++ is 0-based (i=0...n-1)
2. **FFT normalization**: Julia FFTW doesn't normalize by default
3. **Memory management**: Julia garbage-collected vs C++ manual new/delete

### Debugging Tips
- Check FFT normalization: `norm(ifft(fft(x)) - x)` should be ~machine epsilon
- Verify map evaluation: `P(x₀, μ)` should return to Poincare section (q₃=0)
- Monitor small divisors: Plot |k·ω| for all modes to check resonances

## Development Status

**Last Updated**: 2026-08-05

**Current Milestone**: STEP 1 Testing  
**Next Milestone**: STEP 2-4 Implementation

See `docs/JULIA_IMPLEMENTATION_STATUS.md` for detailed progress tracking.

## Contributing

To extend this implementation:

1. **Add STEP 2**: Implement `build_symplectic_frame!()`
   - Fourier derivative: `deriva(paramF, k)`
   - Metric tensor evaluation at each grid point
   - Gram matrix computation

2. **Add STEP 3**: Implement `solve_cohomological_equations!()`
   - Small divisor solver
   - Twist matrix extraction and inversion
   - Average mode handling

3. **Add STEP 4**: Implement `update_parameterization!()`
   - Matrix multiplication at each grid point
   - Fourier conversion

4. **Add continuation loop**: Iterate until convergence
   ```julia
   for iter in 1:maxiter
       error, converged, _ = kam_torus!(paramR, paramF, ω, μ, H)
       if converged
           break
       end
   end
   ```

## Contact

For questions or bug reports:
- See C++ reference: `src/cpp/param.cc`
- Read analysis doc: `docs/analysis/KAM_TORUS_BUG_ANALYSIS.md`
- Check test output: `run_test.jl`
