# Computing Invariant Tori Around Near-Rectilinear Halo Orbits

Computational framework for computing Lagrangian tori (two-dimensional quasi-periodic invariant manifolds) surrounding Near-Rectilinear Halo Orbits (NRHOs) in the Circular Restricted Three-Body Problem (CR3BP), with application to Enceladus mission design.

## Overview

This code implements the **parameterization method** of Haro and Luque for computing invariant tori using:
- **C++** for high-performance numerical computation (KAM torus correction, Poincaré maps, FFT)
- **Julia** for high-level analysis, visualization, and initial condition generation

The Near-Rectilinear Halo Orbits (NRHOs) around Enceladus are a potential mission design solution, and this code computes quasi-periodic orbits (QPOs) around the NRHOs.

## Quick Start

### Build

```bash
make
```

This builds the `param` executable in `bin/`.

### Run

```bash
./bin/param data/initial_conditions/approxQPO.csv
```

This computes a torus correction starting from the initial approximation in `approxQPO.csv`.

### Julia Analysis

```bash
julia --project=. src/julia/main.jl
```

## Directory Structure

```
├── src/                       # Source code
│   ├── cpp/                   # C++ implementation
│   │   ├── param.cc           # Main KAM torus computation
│   │   ├── seccp.c            # Poincaré section integrator
│   │   ├── fluxvp.c           # Flux/flow integrator
│   │   ├── rk78vp.c           # Runge-Kutta integrator
│   │   └── headers/           # C/C++ headers
│   │       ├── complex.h      # Complex number class
│   │       ├── grid.h         # Fourier grid operations
│   │       ├── matrix.h       # Matrix data structures
│   │       └── ...
│   └── julia/                 # Julia analysis code
│       ├── main.jl            # Main interface
│       ├── dynamics.jl        # CR3BP dynamics
│       └── parameterization-method.jl
│
├── data/                      # Input data
│   ├── initial_conditions/    # Initial torus approximations
│   │   ├── approxQPO.csv
│   │   ├── halo.csv
│   │   └── NRHO_*.csv
│   └── config/                # Configuration files
│
├── output/                    # Generated torus output files
├── results/                   # Test results and logs
├── build/                     # Build artifacts (.o files)
├── bin/                       # Executables
├── docs/                      # Documentation
│   ├── analysis/              # Technical analysis
│   │   ├── KAM_TORUS_BUG_ANALYSIS.md
│   │   ├── TEST_RESULTS.md
│   │   └── CONTINUATION_FIX.md
│   ├── paper/                 # SIAM conference paper
│   │   ├── paper.tex
│   │   ├── references.bib
│   │   └── Makefile
│   └── STATUS.md
│
├── figures/                   # Generated figures
├── Makefile                   # Build system
├── Project.toml               # Julia dependencies
└── README.md                  # This file
```

## Dependencies

### C++ Compilation
- GCC/G++ compiler
- Standard C/C++ libraries
- Math library (-lm)

### Julia (Optional for analysis)
- Julia 1.6+
- Packages (see Project.toml):
  - OrdinaryDiffEq
  - StaticArrays
  - DiffEqBase

## Algorithm

The code implements the parameterization method for invariant tori:

1. **Represent torus** as Fourier series: K(θ) where θ are angle variables
2. **Check invariance**: E = Map(K(θ)) - K(θ + ω)
3. **Build symplectic frame**: Tangent (L) and normal (N) bundles
4. **Solve cohomological equations**: Find corrections ξ_L, ξ_N
5. **Update parameterization**: K_new = K + L·ξ_L + N·ξ_N
6. **Iterate** until error < tolerance or max iterations reached

See `docs/paper/paper.tex` for full mathematical description.

## Key Features

- ✅ Adaptive grid refinement (automatically increases resolution if Fourier tail too large)
- ✅ Three methods for normal frame construction (Case 1, 2, 3)
- ✅ Parameter continuation with respect to epsilon
- ✅ Well-conditioned cohomological equation solver
- ✅ Comprehensive error checking and diagnostics

## Bug Fixes (2026-08-05)

Recent improvements:
- Fixed duplicate `kam_torus()` call that caused compounding Fourier truncation errors
- Added continuation loop termination criteria to prevent infinite loops
- Repository reorganized for better structure

See `docs/analysis/` for detailed bug analysis and test results.

## Research Paper

A complete SIAM conference paper is included in `docs/paper/`:

**Title**: "Computing Invariant Tori Around Near-Rectilinear Halo Orbits for Enceladus Mission Design"

Build the paper:
```bash
cd docs/paper
make
```

## References

Based on work by:
- Àlex Haro and Rafael de la Llave (parameterization method)
- Josep-Maria Mondelo (CR3BP applications)

Key paper: Haro, À., et al. (2016). *The parameterization method for invariant manifolds*. Springer.

## License

[Specify license here]

## Contact

Jared Blanchard  
GNC Engineer, True Anomaly  
[Jared.Blanchard@trueanomaly.space](mailto:Jared.Blanchard@trueanomaly.space)

## Citation

If you use this code in your research, please cite:

```bibtex
@misc{blanchard2026nrho,
  author = {Blanchard, Jared},
  title = {Computing Invariant Tori Around Near-Rectilinear Halo Orbits for Enceladus Mission Design},
  year = {2026},
  url = {https://github.com/jared711/nrho-lagrangian-tori}
}
```
