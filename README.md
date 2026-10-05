# Lagrangian tori around Saturn–Enceladus NRHOs

Computes quasi-periodic orbits (QPOs) around stable Near-Rectilinear Halo Orbits (NRHOs) of the
Saturn–Enceladus circular restricted three-body problem, as Lagrangian invariant 2-tori of the
Poincaré map on q3 = 0 at fixed energy, with the parameterization method of Haro and Luque.

**Current status: [`docs/STATUS_2026-10-05.md`](docs/STATUS_2026-10-05.md).**
Who wrote what (Barcelona code, Jared Blanchard's code, AI-assisted changes):
[`docs/PROVENANCE.md`](docs/PROVENANCE.md). History: [`docs/hist/`](docs/hist/README.md).

## Results so far

Five families, 454 tori, around NRHO Ids 23, 37, 53, 85 and 97 of `data/initial_conditions/NRHO_L2.csv`,
with sup-norm invariance errors from ~4e-13 to 1e-11 LU (0.1–2 mm). Torus radii on the section
range from 0.67 to 14.7 km; the perilune altitude of the NRHO limits the useful size (the larger
tori of Id 97, perilune 4.2 km, impact Enceladus). Data in `results/families/`, figures in
`figures/`, paper draft in `docs/paper/`.

## Build and run

Requirements: a C/C++ compiler with OpenMP (GCC; on macOS use Homebrew `gcc` or `libomp`), Julia
(the environment is `Project.toml`/`Manifest.toml` at the repository root), and pdflatex/bibtex
for the paper.

```bash
make                                                   # builds bin/param
./bin/param <anyinput.csv> --test-twist 2              # check: integrable twist map, quadratic convergence
julia --project=. src/julia/family_pipeline.jl 53 /tmp/fam53   # set up a family (writes setup.txt)
```

The full workflow (continuation, analysis, figures) is in the status document.

`bin/param` options: `--domega d1 d2`, `--eps-max`, `--max-steps`, `--tol-floor` (continuation),
`--free-omega`, `--no-filter`, `--no-real`, `--fixed-frame`, `--tol-coho` (method variants),
`--test-twist`, `--fdcheck`, `--orbit`, `--eval` (checks and helpers). Threads: `OMP_NUM_THREADS`.

## Layout

```
src/cpp/            param.cc (kam_torus, map_CR3BP, continuation), RTBP integrator and Poincaré
                    section library (rtbphp, seccp, fluxvp, rk78vp, ...), headers/ (grid/FFT, matrix)
src/julia/          analysis and set-up tools (see src/julia/README.md)
data/               NRHO initial conditions (NRHO_L1.csv, NRHO_L2.csv), older inputs
results/families/   computed families (one directory per NRHO; READMEs give commands and code versions)
figures/            family and QPO figures
docs/               STATUS_<date>.md, PROVENANCE.md, paper/, hist/
```

## Method

For a torus K: T² → section with rotation vector ω, solve F(K(θ)) = K(θ + ω) by a quasi-Newton
method: tangent frame L = DK, symplectic normal frame N, torsion T, and two cohomological
equations solved in Fourier space (Figueras, Haro & Luque 2017; Haro et al. 2016, ch. 4). The
CR3BP-specific adaptations (torus-adapted symplectic coordinates, sup-norm error, filtering,
anisotropic refinement, continuation in ω) are described in the paper and in the commit messages.

## Credits

The KAM torus code is by Àlex Haro and Alejandro Luque (2014); the RTBP and Poincaré-section
library by Àlex Haro and Josep-Maria Mondelo (2016, LGPL). The CR3BP adaptation and the Julia code
are by Jared Blanchard, with AI-assisted changes in 2026 recorded in `docs/PROVENANCE.md`.

## Contact

Jared Blanchard — [Jared.Blanchard@trueanomaly.space](mailto:Jared.Blanchard@trueanomaly.space)
