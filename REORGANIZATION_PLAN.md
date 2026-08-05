# Repository Reorganization Plan

## Current State (Messy)
- All source files mixed in root directory
- C/C++ files, Julia files, headers, data files all together
- Documentation scattered
- Test outputs in root
- Object files (.o) in root

## Proposed Structure

```
nrho-lagrangian-tori/
├── README.md                      # Main project README
├── LICENSE                        # License file (if applicable)
├── .gitignore                     # Git ignore (updated)
│
├── src/                          # Source code
│   ├── cpp/                      # C++ implementation
│   │   ├── param.cc              # Main KAM torus computation
│   │   ├── seccp.c               # Poincaré section integrator
│   │   ├── fluxvp.c              # Flux/flow integrator
│   │   ├── rk78vp.c              # Runge-Kutta integrator
│   │   ├── rtbphp.c              # Hamiltonian computation
│   │   ├── campvp.c              # CR3BP vector field
│   │   ├── scread.c              # Input file reader
│   │   ├── vbprintf.c            # Verbose printing utilities
│   │   └── headers/              # C/C++ headers
│   │       ├── complex.h
│   │       ├── grid.h            # Fourier grid operations
│   │       ├── matrix.h          # Matrix data structures
│   │       ├── utils.h
│   │       ├── seccp.h
│   │       ├── fluxvp.h
│   │       ├── rk78vp.h
│   │       ├── rtbphp.h
│   │       ├── campvp.h
│   │       ├── scread.h
│   │       └── vbprintf.h
│   │
│   └── julia/                    # Julia analysis code
│       ├── main.jl               # Main Julia interface
│       ├── dynamics.jl           # CR3BP dynamics
│       ├── parameterization-method.jl  # Torus utilities
│       └── util.jl               # Utility functions
│
├── data/                         # Input data files
│   ├── initial_conditions/       # IC files
│   │   ├── approxQPO.csv         # Initial torus approximation
│   │   ├── approxQPO_p128.csv    # High-res version
│   │   ├── halo.csv              # Halo orbit
│   │   ├── NRHO_L1.csv           # L1 NRHO
│   │   └── NRHO_L2.csv           # L2 NRHO
│   │
│   └── config/                   # Configuration files
│       ├── curve_start_names.txt
│       ├── nrho_start_names.txt
│       ├── nrho_start.txt
│       └── torus_start_names.txt
│
├── output/                       # Output files (gitignored)
│   ├── .gitkeep                  # Keep directory in git
│   └── README.md                 # Explains output file formats
│
├── results/                      # Test results and logs
│   ├── param_single_kam_test.txt
│   ├── param_output.txt
│   ├── full_param_output.txt
│   ├── julia_output.txt
│   ├── julia_run2.txt
│   ├── julia_run3.txt
│   └── mytraj.txt
│
├── build/                        # Build artifacts (gitignored)
│   ├── .gitkeep
│   └── *.o files go here
│
├── bin/                          # Executables (gitignored)
│   ├── .gitkeep
│   └── param (executable goes here)
│
├── docs/                         # Documentation
│   ├── analysis/                 # Technical analysis docs
│   │   ├── KAM_TORUS_BUG_ANALYSIS.md
│   │   ├── TEST_RESULTS.md
│   │   ├── CONTINUATION_FIX.md
│   │   └── SESSION_SUMMARY_2026-08-05.md
│   │
│   ├── paper/                    # Conference paper
│   │   ├── paper.tex
│   │   ├── references.bib
│   │   ├── Makefile
│   │   ├── README.md
│   │   ├── PAPER_SUMMARY.md
│   │   └── figures/
│   │       └── README.md
│   │
│   └── STATUS.md                 # Project status
│
├── figures/                      # Generated figures
│   └── 2torusIC.svg
│
├── scripts/                      # Utility scripts
│   └── README.md                 # For future build/test scripts
│
├── tests/                        # Test files
│   ├── test.cpp
│   └── README.md                 # For future test suite
│
├── Makefile                      # Updated build system
├── Project.toml                  # Julia project file (if exists)
└── Manifest.toml                 # Julia manifest
```

## Benefits

1. **Clear Separation**: Source code separated by language (C++/Julia)
2. **Clean Root**: Only essential files in root (README, Makefile, license)
3. **Build Artifacts Isolated**: .o files and executables in separate dirs
4. **Data Organization**: Input data separated from results/output
5. **Documentation Grouped**: All docs in one place with subcategories
6. **Gitignore Friendly**: Build/output dirs can be easily ignored
7. **Professional Structure**: Follows standard open-source conventions
8. **Scalable**: Easy to add tests, benchmarks, examples later

## Migration Steps

1. Create directory structure
2. Move files to new locations
3. Update Makefile for new paths
4. Update .gitignore
5. Update Julia imports/includes
6. Update README with new structure
7. Test build system
8. Test Julia code
9. Commit reorganization

## Files to Keep in Root

- `README.md` - Main project description
- `Makefile` - Build system
- `LICENSE` - License file
- `.gitignore` - Git configuration
- `Project.toml` / `Manifest.toml` - Julia dependencies
- `.gdbinit` - Debugger config (optional, could move to scripts/)

## Files to Delete/Clean

- `*.o` files - Will be regenerated in build/
- `mypoint.txt` - Temporary test file (empty)
- Old test outputs - Move to results/ or delete if obsolete
- `param` executable - Will be regenerated in bin/
