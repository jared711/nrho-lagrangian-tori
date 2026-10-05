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

## How to run everything

All commands are run from the repository root unless a `cd` is shown. `REPO` stands for the
absolute path of the repository (`REPO=$(pwd)`).

### 0. Requirements

* C/C++ compiler with OpenMP: GCC on Linux; on macOS install Homebrew `gcc` and build with
  `make CC=gcc-14 CXX=g++-14` (or the installed version).
* Julia (tested with 1.12). The environment is `Project.toml`/`Manifest.toml` at the repository
  root; the first time, run `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
* pdflatex and bibtex for the paper.

### 1. Build and check

```bash
make                                                     # builds bin/param (-O2, OpenMP)
./bin/param data/initial_conditions/approxQPO.csv --test-twist 2   # twist map: quadratic convergence to ~1e-16
./bin/param data/initial_conditions/approxQPO.csv --fdcheck        # Jacobian vs finite differences (~6e-7)
```

`OMP_NUM_THREADS=<n>` sets the number of threads used for the map evaluations.

### 2. Set up a family around an NRHO

Pick an elliptic member of `data/initial_conditions/NRHO_L2.csv` (see
`results/families/elliptic_members_scan.txt` for perilune altitudes and resonance distances), e.g.
Id 53:

```bash
julia --project=. src/julia/family_pipeline.jl 53 /tmp/fam53
cat /tmp/fam53/setup.txt        # rotation numbers, domega, and the param commands
```

This writes `/tmp/fam53/start.csv` (a 32×32 start torus fitted to an orbit of the NRHO neighbourhood) and
`results/families/nrho_id53/nrho_fixed_point.txt` (needed later by the visualization).

### 3. Continue the family (one directory per branch)

`bin/param` writes each converged torus as `output_torus<ε>` into the directory it runs in. The
continuation step is column 11 of the first line of the start file: positive grows the tori
(outward), negative shrinks them toward the NRHO (inward).

```bash
D1D2="2.8729e-4 4.2791e-4"      # the "domega" line of setup.txt
mkdir -p /tmp/fam53_out /tmp/fam53_in
awk 'NR==1{$11="0.01"}  {print}' /tmp/fam53/start.csv > /tmp/fam53_out/start.csv
awk 'NR==1{$11="-0.01"} {print}' /tmp/fam53/start.csv > /tmp/fam53_in/start.csv
cd /tmp/fam53_out && OMP_NUM_THREADS=12 nohup $REPO/bin/param start.csv --domega $D1D2 \
     --eps-max 100  --max-steps 3000 --tol-floor 1e-11 > fam.log 2>&1 &
cd /tmp/fam53_in  && OMP_NUM_THREADS=12 nohup $REPO/bin/param start.csv --domega $D1D2 \
     --eps-max -0.95 --max-steps 3000 --tol-floor 1e-11 > fam.log 2>&1 &
```

* Progress: `grep -E '^-?[0-9]' fam.log` lists the accepted tori (ε, grid, …, error).
* A branch stops when the continuation step falls below its minimum, at `--eps-max`, or after
  `--max-steps` steps. The grid is refined automatically (only in the angle that needs it) when
  Newton stalls near the tolerance.
* `--tol-floor` is the accepted sup-norm invariance error in LU (1e-11 LU ≈ 2.4 mm; relax it, e.g.
  to 1e-10, to reach larger tori).
* To extend a branch from its last torus (ε is then relative to that torus):

  ```bash
  julia --project=. src/julia/restart_from.jl /tmp/fam53_out/output_torus<last ε> /tmp/fam53_out2/start.csv 0.01
  ```

### 4. Save a family

Copy each branch into `results/families/nrho_id<Id>/<branch>/` (tori, `start.csv`, the log as
`run.log`), gzip large files (`gzip -9 output_torus* run.log`; the Julia tools read `.gz`
directly), and add a README with the command and the code version.

### 5. Tables and figures

```bash
B=results/families/nrho_id53
NRHO_ID=53 RHO_LINEAR=<rho_lin_1>,<rho_lin_2> julia --project=. src/julia/family_analysis.jl \
     figures/family_id53/id53 inward=$B/inward_1,$B/inward_2 outward=$B/outward_1,$B/outward_2
NRHO_ID=53 julia --project=. src/julia/qpo_visualization.jl figures/qpo3d/id53 30 \
     <smallest torus file> <middle torus file> <largest torus file>
```

`family_analysis.jl` writes a CSV table (ε, rotation numbers, radius, grid, error per torus) and
four figures; several directories of one branch are joined with commas. `qpo_visualization.jl`
integrates the QPOs, prints their distance to the NRHO and minimum altitude above Enceladus, and
writes a 3D view, the distance to the NRHO over time, and apolune cross-sections (`NCROSS`
revolutions, default 300).

### 6. Paper

```bash
cd docs/paper && make           # paper.pdf
```

### `bin/param` reference

```
bin/param <input.csv> [options]
  continuation:  --domega d1 d2   --eps-max E   --max-steps N   --tol-floor T (LU)
  method:        --free-omega  --no-filter  --no-real  --fixed-frame  --tol-coho t
  checks/tools:  --test-twist [Case]  --fdcheck  --orbit niter q1 q2 p1 p2  --eval
```

Input file: first line `toltail tolinva tolinte ω1 ω2 ε H μ n1 n2 Δε 1`, then one line per grid
point `i1 i2 q1 q2 p1 p2` (1-based indices, second index fastest). Output files have the same
values, one per line in the header and 0-based indices.

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
