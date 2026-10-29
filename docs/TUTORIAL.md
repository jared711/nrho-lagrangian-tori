# Tutorial: computing Lagrangian tori around Saturn–Enceladus NRHOs

This tutorial walks through the code in this repository from a fresh clone to a family of tori
and its figures. The [README](../README.md) is the short command reference; this document
explains what each step does, what the output means, and what to do when something goes wrong.

All commands run from the repository root unless a `cd` is shown. `REPO=$(pwd)` is the absolute
path of the repository. Output shown below was produced with the code at commit `932e746` on a
2-core Linux machine.

Contents

1. [What the code computes](#1-what-the-code-computes)
2. [Build and check](#2-build-and-check)
3. [Conventions: coordinates, units, files](#3-conventions-coordinates-units-files)
4. [Walkthrough A: continue a family from a saved start torus (no Julia needed)](#4-walkthrough-a-continue-a-family-from-a-saved-start-torus)
5. [Reading the log](#5-reading-the-log)
6. [Walkthrough B: a new family around any elliptic NRHO](#6-walkthrough-b-a-new-family-around-any-elliptic-nrho)
7. [Production runs: both branches, restarts, saving](#7-production-runs-both-branches-restarts-saving)
8. [Tables, figures and QPOs](#8-tables-figures-and-qpos)
9. [How large can the tori get? (regular-region scans)](#9-how-large-can-the-tori-get)
10. [`bin/param` reference](#10-binparam-reference)
11. [Troubleshooting](#11-troubleshooting)
12. [Map of the source code](#12-map-of-the-source-code)
13. [Building the paper](#13-building-the-paper)

---

## 1. What the code computes

An NRHO that is elliptic in both normal directions is surrounded by quasi-periodic orbits (QPOs)
lying on invariant 3-tori of the flow. Cutting the flow with the Poincaré section q3 = 0 (upward
crossings) at fixed energy H turns the NRHO into a fixed point z\* of a 4-dimensional symplectic
map F, and each 3-torus of the flow into an invariant **2-torus** of F. The code finds such a
2-torus as a parameterization K: T² → R⁴ that satisfies the invariance equation

    F(K(θ)) = K(θ + ω),      θ ∈ T² = [0,1)²,

where ω = (ω₁, ω₂) is the rotation vector of the torus (in cycles per map iterate). K is stored
as its values on an n₁ × n₂ grid of angles; Fourier transforms move between grid values and
Fourier coefficients.

The pipeline:

```
NRHO_L2.csv ──family_pipeline.jl──▶ start.csv (32×32 torus fitted to an orbit) + setup.txt
                                       │
                                       ▼
                      bin/param  (Newton on the invariance equation,
                                  continuation in ω, grid refinement)
                                       │
                                       ▼
                      output_torus<ε> files, one per converged torus
                                       │
             ┌─────────────────────────┴───────────────────────┐
             ▼                                                 ▼
   family_analysis.jl (table, figures)            qpo_visualization.jl (integrate
                                                   QPOs, altitude above Enceladus)
```

The heavy lifting is in C++ (`bin/param`); the set-up and analysis are in Julia. The solver needs
no Julia once you have a start file, so you can run continuations on any machine with a C++
compiler.

## 2. Build and check

### Requirements

* A C/C++ compiler with OpenMP. On Linux, GCC. On macOS, Apple clang has no OpenMP: install
  Homebrew `gcc` and build with `make CC=gcc-14 CXX=g++-14` (or the version you have).
* Julia (tested with 1.12; `Manifest.toml` pins 1.12.5) for set-up and analysis. The first time:
  `julia --project=. -e 'using Pkg; Pkg.instantiate()'`. `ThreeBodyProblem.jl` is one of the
  dependencies.
* pdflatex and bibtex for the paper.

### Build

```bash
make            # builds bin/param with -O2 and OpenMP
```

A few `-Wmaybe-uninitialized` warnings from `param.cc` are expected and harmless.

### Check 1: the twist map (does the Newton method work?)

```bash
./bin/param data/initial_conditions/approxQPO.csv --test-twist 2 | grep -E "Newton iteration|Error of invariance"
```

This runs the torus solver on an integrable twist map whose invariant tori are known exactly,
starting from a deliberately wrong torus (1% wrong action plus a 1e-4 harmonic). It uses only
the grid size and tolerances of the input file, not the CR3BP. Expected output:

```
# Newton iteration 1
#     - Error of invariance: 2.124e-04 (max over grid, physical), ...
# Newton iteration 2
#     - Error of invariance: 1.269e-06 (max over grid, physical), ...
# Newton iteration 3
#     - Error of invariance: 7.137e-11 (max over grid, physical), ...
# Newton iteration 4
#     - Error of invariance: 1.249e-16 (max over grid, physical), ...
```

The error roughly squares at each step (2e-4 → 1e-6 → 7e-11 → 1e-16): quadratic convergence. If
you change anything in `kam_torus`, run this first; losing quadratic convergence here means the
change broke the method, not the CR3BP set-up.

### Check 2: the differential of the Poincaré map

```bash
./bin/param data/initial_conditions/approxQPO.csv --fdcheck | tail -2
```

```
# max |analytic - fd| / max(1,|fd|) = 5.831660621881545e-07
# map failures: 0
```

The analytic DF (from the variational equations) agrees with central finite differences to about
6e-7, which is the accuracy of the finite differences themselves. A value near 1 means the
Jacobian and the map disagree, and Newton will not converge on the CR3BP.

### Threads

The N map evaluations of each Newton step run in parallel with OpenMP. Set
`OMP_NUM_THREADS=<n>` (the number of physical cores is a good choice).

## 3. Conventions: coordinates, units, files

### Model and units

* Saturn–Enceladus CR3BP, μ = 1.901109735892602e-7.
* Length unit LU = 238 042 km (the Saturn–Enceladus distance). So 1e-11 LU ≈ 2.4 mm, and
  1e-5 LU ≈ 2.4 km. Enceladus' mean radius is 252.1 km.
* **Barcelona convention.** The solver uses the Hamiltonian form of the CR3BP of the Barcelona
  RTBP library (`src/cpp/rtbphp.*`): positions q, momenta p, and

      H = ½|p|² + q₂p₁ − q₁p₂ − (1−μ)/r₁ − μ/r₂,

  with Saturn at q₁ = μ and Enceladus at q₁ = μ − 1. That is the standard rotating frame
  rotated by 180° about the z-axis, so Enceladus sits near q₁ ≈ −1. Energy and Jacobi constant
  are related by H = −C/2. `src/julia/util.jl` converts between the two: `rv2pq` (standard
  position and velocity → Barcelona (q, p)) and `pq2rv`.
* **Section:** q₃ = 0, upward crossing. A point on the section is z = (q₁, q₂, p₁, p₂); p₃ is
  recovered from H. Every torus file, start file and fixed point uses these four coordinates.

### Input file (`start.csv`)

The first line is a header with 12 numbers. After it comes one line per grid point.

| field | name | meaning |
|---|---|---|
| 1 | `toltail` | tail tolerance in LU: if the Fourier coefficients in the upper part of the spectrum exceed it, the grid is refined (normally 1e-11) |
| 2 | `tolinva` | target invariance error for Newton (normally 1e-12; in practice `--tol-floor` decides acceptance, see §5) |
| 3 | `tolinte` | tolerance for intermediate linear algebra (1e-17) |
| 4–5 | ω₁ ω₂ | rotation vector of the torus (cycles per iterate) |
| 6 | ε | continuation parameter of this torus (usually 0) |
| 7 | H | energy, Barcelona convention |
| 8 | μ | mass parameter |
| 9–10 | n₁ n₂ | grid size |
| 11 | Δε | **first continuation step: its sign sets the direction** (+ outward, − inward) |
| 12 | `1` | flag, always 1 |

Each grid point line is `i1 i2 q1 q2 p1 p2`. The indices are 1-based and the second index runs
fastest.

### Output files (`output_torus<ε>`)

`bin/param` writes one file per converged torus into the **current directory**. The filename
holds ε with 6 decimals, for example `output_torus0.020000`. The file carries the same
information as an input file, laid out differently:

* the 12 header values one per line instead of on one line;
* 0-based indices;
* field 11 holds the step size the continuation would use next.

**Caveat:** field 1 of an output file is the tail tolerance after it has been rescaled to the
solver's local coordinates (for example 8.4e-08), not the 1e-11 LU you started with. Reset it
when you turn an output file into a new start file (see §7.3).

### Other files

* `data/initial_conditions/NRHO_L2.csv`: the L2 halo/NRHO family. Each row has an Id, the
  initial state in the standard rotating frame, the Jacobi constant, the period and the
  stability index.
* `results/families/nrho_id<Id>/nrho_fixed_point.txt`: one line, `q1 q2 p1 p2 H μ T`. These
  are the fixed point z\* of the section map, the energy, μ, and the NRHO period in TU. It is
  written by `family_pipeline.jl` and read by the visualization.
* `results/families/elliptic_members_scan.txt`: for every elliptic NRHO member, the perilune
  altitude above Enceladus, the linear rotation numbers ρ, and the nearest low-order resonance.
  Use it to choose a member (§6).

## 4. Walkthrough A: continue a family from a saved start torus

This walkthrough needs only `bin/param`. It reproduces the first steps of the published Id 53
family.

```bash
REPO=$(pwd)
mkdir -p /tmp/demo53 && cd /tmp/demo53
cp $REPO/results/families/nrho_id53/outward_1/start.csv .
head -1 start.csv
```

```
9.9999999999999994e-12 9.9999999999999998e-13 1.0000000000000001e-17 0.33938306646592242 0.15261138018692144 0.0 -1.5000173895290601 1.9011097358926021e-07 32 32 0.01 1
```

This header describes a 32×32 torus with ω = (0.339383…, 0.152611…) at H = −1.50001738952906,
with a first step of +0.01 (outward). Run five continuation steps:

```bash
OMP_NUM_THREADS=2 $REPO/bin/param start.csv --domega 2.8729e-4 4.2791e-4 \
    --eps-max 100 --max-steps 5 --tol-floor 1e-11 > fam.log 2>&1
grep -E '^-?[0-9]' fam.log
ls
```

```
0.000000000000000e+00 32 32 2.030196607732267e+04 8.766686569272552e-01 6.796245161319564e-01 5.821053850716433e-13
2.000000000000000e-02 32 32 2.034855099885603e+04 8.861395962042143e-01 6.897589285487148e-01 8.968713227939802e-13
4.666666666666666e-02 32 32 1.999903479678307e+04 8.894072819046436e-01 6.959823202017972e-01 7.536081889551401e-13
9.999999999999999e-02 32 32 1.957208461789557e+04 9.015972286288051e-01 7.129269453793711e-01 5.215200012849087e-13
1.711111111111111e-01 32 32 1.872674048357088e+04 9.091202285606091e-01 7.286427697455350e-01 6.000615339978986e-13

fam.log  output_torus0.000000  output_torus0.020000  output_torus0.046667  output_torus0.100000  output_torus0.171111  start.csv
```

The run took 6 seconds on 2 cores. Here is what happened:

1. **Step 1 (ε = 0)** re-converges the start torus itself. Its error drops from about 8e-12
   (orbit fit) to 5.8e-13 LU.
2. **Steps 2–5** move along the family. The rotation vector follows

       ω(ε) = ω_start + ε · Δω,     Δω = (2.8729e-4, 4.2791e-4)  (from --domega)

   Δω is chosen so that ε = −1 lands on the linear rotation numbers of the NRHO itself: ε < 0
   shrinks the tori toward the NRHO, and ε > 0 grows them.
3. **Step size control.** Each new torus is predicted by a secant through the last two converged
   tori and then corrected by Newton. After a success the step is multiplied by 4 / (number of
   Newton iterations), clamped to [½, 2]. That is why the steps grew 0.02 → 0.027 → 0.053 →
   0.071. After a failure the step is halved.

The columns of each accepted-torus line:

| col | meaning |
|---|---|
| 1 | ε |
| 2–3 | grid n₁ n₂ |
| 4 | sum of the entries of the inverse average torsion ⟨T⟩⁻¹, in the solver's local coordinates (printed as `Norm inverse twist` in the Newton lines). Its magnitude grows when the twist condition weakens |
| 5–6 | Sobolev-type norms of the q₁ and q₂ components of K in local coordinates (they track the torus size and smoothness, not physical km) |
| 7 | **sup-norm invariance error max\|F(K(θ)) − K(θ+ω)\| over the grid, in LU**. This is the number that matters |

The torus radius in km is computed afterward by `family_analysis.jl` (§8).

## 5. Reading the log

`fam.log` is verbose: lines starting with `#` describe each Newton step, and lines starting with a
number are accepted tori. Useful filters:

```bash
grep -E '^-?[0-9]' fam.log                     # accepted tori
grep -E 'Error of invariance' fam.log | tail   # Newton progress
grep -E '^# [A-Z]' fam.log | grep -v Newton\ iteration | sort | uniq -c   # events
```

What one Newton step reports (`# Newton iteration k`):

* `Error of invariance: a (max over grid, physical), b (Fourier l1, local)`: a is the sup-norm
  error in LU, the one used for acceptance.
* `Error spectrum near the filter cutoff: x% (theta1), y% (theta2)`: how much of the error sits
  in the highest retained Fourier modes. Above ~20% in an angle means the grid is too coarse in
  that angle.
* `Norm inverse twist`: the sum of the entries of ⟨T⟩⁻¹ (signed, so it can be negative). The twist (non-degeneracy) condition fails when its magnitude blows up.
* `Lagrangian defect |L^T Omega L|`: should be small. It measures how far the torus is from
  Lagrangian, which is what the approximately symplectic frame relies on.
* `Frame check |L^T Omega N + I|`: should be at round-off (~1e-15).
* `max |Im K| after the step: … (removed)`: K is real, so the round-off imaginary part is
  dropped (`--no-real` keeps it).

Events and what they mean:

| message | meaning | what happens next |
|---|---|---|
| `Newton stalled; best error e (physical), accepted (tol-floor T)` | Newton stopped improving (a step reduced the error by less than 10%) at the noise floor of the map, and the best error is below `--tol-floor` | the best torus is accepted |
| `… not accepted` | stalled above `--tol-floor` | the step is halved |
| `Newton stalled above tol-floor; refining the grid to n₁ x n₂` | stalled within 100× of the floor, and the error spectrum shows truncation | the grid is doubled in the truncated angle(s) only, and Newton restarts |
| `Accepted torus close to tol-floor and truncated; refining the grid` | accepted with error > 0.3 · tol-floor and truncated | proactive refinement before the next step, so the next step has margin |
| `We need more Fourier modes` | the tails exceed `toltail` | the grid is doubled in the angle(s) with large tails |
| `We reched the maximum number of Fourier modes` | n₁·n₂ would exceed 262 144 (512×512) | the run ends |
| `Newton method does not converge; halving the continuation step` | Newton failed (10 iterations or diverged) | step halved; **the run ends silently once the step is below 1/1000 of the initial step** |
| `Reached maximum epsilon` / `Reached maximum continuation steps` | `--eps-max` or `--max-steps` reached | the run ends |

**Why the stall detection exists.** Each map evaluation integrates the flow numerically, so F
itself carries an error of roughly 1e-13 to 1e-12 LU (the "noise floor"). Once the torus error
reaches that floor, further Newton steps only amplify the noise through the small divisors of
the cohomological equations. The code therefore keeps the best torus seen and stops as soon as
improvement stalls. **`--tol-floor` is the real acceptance threshold**, and `tolinva` in the
header is usually never reached.

## 6. Walkthrough B: a new family around any elliptic NRHO

### 6.1 Choose a member

```bash
head -20 results/families/elliptic_members_scan.txt
```

```
Id 111  T=2.2037  perilune altitude     0.2 km   rho=(0.362143, 0.202330) nearest res k=(0, -5) dist=0.0117
Id 109  T=2.2053  perilune altitude     0.7 km   rho=(0.360906, 0.201352) nearest res k=(-5, -1) dist=0.0059
Id 101  T=2.2121  perilune altitude     3.1 km   rho=(0.356035, 0.196740) nearest res k=(0, -5) dist=0.0163
...
```

There are two things to look at:

* **Perilune altitude.** QPOs swing around the NRHO, so their closest approach to Enceladus is
  lower than the NRHO's perilune. Around Id 97 (4.2 km) the larger tori hit the surface, while
  around Id 53 (17.1 km) the largest torus still clears it by 8.9 km.
* **Distance to the nearest resonance** k·ρ ∈ Z (low order |k|). Close to a low-order resonance
  the tori break up early and the cohomological equations have small divisors.

### 6.2 Run the set-up

```bash
julia --project=. src/julia/family_pipeline.jl 53 /tmp/fam53
```

What it does (all via `bin/param --orbit` and `--eval`, driven from `src/julia/pmap_tools.jl`):

1. Reads row Id of `NRHO_L2.csv` and converts it to Barcelona coordinates. H = −C/2.
2. Newton-solves for the **fixed point** z\* of the section map (the NRHO). It checks that all
   four eigenvalues of DF(z\*) lie on the unit circle (elliptic–elliptic) and stops otherwise.
3. Takes the **linear rotation numbers** ρ_lin and the elliptic eigenvectors V₁, V₂ from DF(z\*).
4. Starts an orbit at z\* + a(c₁ Re V₁ + c₂ Re V₂) for a = 4e-5, 2e-5, 1e-5 LU and iterates the map
   8192 times. It accepts the first orbit that is **regular**: the rotation numbers from NAFF
   frequency analysis of the two halves of the orbit agree to 1e-8.
5. **Fits a torus to the orbit.** The orbit fills a 2-torus densely, so the code fits a Fourier
   series (|mᵢ| ≤ 10) to 6000 iterates using their known angles θ = kρ. It requires a fit residual
   below 1e-9 and evaluates the fit on a 32×32 grid.
6. Sets **Δω = ρ − ρ_lin**, so that ε = −1 corresponds to the NRHO.

It writes:

* `/tmp/fam53/start.csv`: the start torus, in input format.
* `/tmp/fam53/setup.txt`: Id, H, μ, period, ρ_lin, ρ, the `domega` line, and the two `param`
  commands.
* `results/families/nrho_id53/nrho_fixed_point.txt`.

The start file's header has Δε = 5e-5 by default. Set the direction and size of the first step
before running (§7.1).

**Mode weights.** `family_pipeline.jl 51 /tmp/fam51 0.3 1` starts the torus with relative
amplitudes (c₁, c₂) = (0.3, 1) in the two elliptic modes. Around some members the regular region
reaches much farther in one mode than in the other (§9), so a different mix gives larger tori.

If it fails with `no regular start torus found`, the member is too close to a resonance or to
the edge of stability at those amplitudes. Try another member or another mode mix.

## 7. Production runs: both branches, restarts, saving

### 7.1 Two branches, one directory each

`bin/param` writes its tori into the directory it runs in, so give each branch its own:

```bash
REPO=$(pwd)
D1D2="2.8729e-4 4.2791e-4"            # the domega line of setup.txt
mkdir -p /tmp/fam53_out /tmp/fam53_in
awk 'NR==1{$11="0.01"}  {print}' /tmp/fam53/start.csv > /tmp/fam53_out/start.csv
awk 'NR==1{$11="-0.01"} {print}' /tmp/fam53/start.csv > /tmp/fam53_in/start.csv

cd /tmp/fam53_out && OMP_NUM_THREADS=12 nohup $REPO/bin/param start.csv --domega $D1D2 \
     --eps-max 100   --max-steps 3000 --tol-floor 1e-11 > fam.log 2>&1 &
cd /tmp/fam53_in  && OMP_NUM_THREADS=12 nohup $REPO/bin/param start.csv --domega $D1D2 \
     --eps-max -0.95 --max-steps 3000 --tol-floor 1e-11 > fam.log 2>&1 &
```

* **Always pass `--eps-max`, `--max-steps` and `--tol-floor`.** Their defaults (0.01, 10 and
  1e-10) are set for testing, so a run without them stops after a few tori.
* Inward, `--eps-max -0.95` stops just short of the NRHO (ε = −1), where the torus degenerates
  to a point.
* Monitor progress with `grep -E '^-?[0-9]' fam.log | tail`.
* Expect the grid to grow (32 → 64 → 128) as the tori grow. Each doubling costs about 4× per
  step: on a 128×128 grid a step takes minutes, not seconds.

### 7.2 When does a branch end?

There are four reasons: `--eps-max`, `--max-steps`, the maximum grid, or the step falling below
1/1000 of its initial size after repeated failures. The last one is the usual end of an outward
branch, and it does not print a final message. The log just ends after a run of `halving the
continuation step` lines. Near the end of a family this is expected and physical. §9 and the
paper's Sec. 6 show that the families end close to the breakdown of the tori.

### 7.3 Restarting or extending a branch from its last torus

```bash
julia --project=. src/julia/restart_from.jl /tmp/fam53_out/output_torus3.271000 /tmp/fam53_out2/start.csv 0.01
```

`restart_from.jl` converts an output file to an input file and resets ε to 0, so ε in the new
run is **relative to this torus**. Keep track of the cumulative ε yourself; see the table in
`results/families/nrho_id53/README.md`. The ω of the restarted torus is its current ω, and Δω
stays the same.

**Reset the tail tolerance.** As noted in §3, `restart_from.jl` copies field 1 unchanged, so the
new start file carries the rescaled tail tolerance (~1e-7) instead of 1e-11. Fix it:

```bash
awk 'NR==1{$1="1e-11"} {print}' start.csv > tmp && mv tmp start.csv
```

**Without Julia**, this awk command does the whole conversion, including the reset:

```bash
awk -v step=0.01 'NR<=12{h[NR]=$1; if(NR==12){h[1]="1e-11"; h[6]="0.0"; h[11]=step;
    s=h[1]; for(i=2;i<=12;i++) s=s" "h[i]; print s}; next}
    {print $1+1, $2+1, $3, $4, $5, $6}' output_torus0.171111 > start.csv
```

(Tested: restarting the walkthrough run from `output_torus0.171111` re-converges it to 5.9e-13 LU
and continues.)

Two reasons to restart:

* **Relax the tolerance to go further.** With `--tol-floor 1e-10` (2.4 cm) the Id 53 outward
  branch reached 15.4 km instead of 14.7 km (`outward_ext*`).
* **Start over on a coarser grid** after a branch refined itself to 128×128 early: resample, then
  continue.

### 7.4 Saving a family

Copy each branch into `results/families/nrho_id<Id>/<branch>/` (`inward_1`, `outward_1`, …;
pieces of one branch get increasing numbers). Each directory gets:

* the `output_torus*` files;
* `start.csv`;
* the log, renamed `run.log`.

Then compress, since the Julia tools read `.gz` directly:

```bash
gzip -9 output_torus* run.log
```

Add a `README.md` to `nrho_id<Id>/` with the member, H, ρ_lin, Δω, the exact command, the code
commit, and the cumulative ε of each piece. The existing READMEs are templates.

## 8. Tables, figures and QPOs

### 8.1 Family table and figures

```bash
B=results/families/nrho_id53
NRHO_ID=53 RHO_LINEAR=0.33909577956895687,0.1521834704601715 \
  julia --project=. src/julia/family_analysis.jl figures/family_id53/id53 \
  inward=$B/inward_1,$B/inward_2,$B/inward_3 outward=$B/outward_1,$B/outward_2,$B/outward_3
```

* Arguments are `<label>=<dir>[,<dir>…]`, one label per branch. Several pieces of one branch are
  joined with commas.
* `RHO_LINEAR` is ρ_lin from `setup.txt` or the family README, and marks the NRHO on the frequency
  plot.
* It writes `<prefix>_table.csv` with columns `branch, eps, omega1, omega2, radius_km, grid,
  error_LU`, and figures `<prefix>_*.png`: rotation numbers, radius vs ρ₂, error vs ρ₂, and
  section views of selected tori.
* It prints a per-branch summary: number of tori, ρ₂ range, radius range, and median and maximum
  error.
* **Radius** is the maximum distance of the grid points from the torus centre in the (q₁, q₂)
  plane of the section, converted to km. This is the "radius on the section" quoted in the paper.

### 8.2 QPOs in configuration space

```bash
NRHO_ID=53 julia --project=. src/julia/qpo_visualization.jl figures/qpo3d/id53 30 \
  $B/inward_3/output_torus-0.071xxx.gz $B/outward_1/output_torus1.2xxxxx.gz $B/outward_3/output_torus0.21xxxx.gz
```

Pick the actual filenames with `ls`. For each torus, the script lifts grid points to the full
6-D state (p₃ from H), integrates the CR3BP for `<nrevs>` NRHO periods, and converts to
Enceladus-centred km in the standard rotating frame. It then:

* prints the **maximum distance from the NRHO** and the **minimum altitude above Enceladus**
  (radius 252.1 km). The second is the number that decides whether a torus is usable for an
  orbiter;
* writes a 3-D view, distance-to-NRHO vs time, and apolune cross-sections (`NCROSS` revolutions,
  default 300).

The NRHO is read from `results/families/nrho_id$NRHO_ID/nrho_fixed_point.txt` (the default
NRHO_ID is 97). `NRHO_FIXED=<file>` points to another fixed-point file.

## 9. How large can the tori get?

These two scripts produced the paper's section on the limits of the families:

```bash
julia --project=. src/julia/stability_edge.jl 53 51 97     # regular-region extent per member
julia --project=. src/julia/breakdown_analysis.jl figures/breakdown/limits
```

`stability_edge.jl` launches orbits at z\* + a(c₁ Re V₁ + c₂ Re V₂) along three mode mixes
(1, 1), (1, 0.3) and (0.3, 1), with growing a. An orbit counts as regular while the NAFF rotation
numbers of its two halves (4096 iterates) agree to 1e-7. For each direction it reports the
largest regular amplitude, that orbit's maximum distance from the NRHO on the section, and the
first irregular amplitude. Results are in `results/families/stability_edge_scan.txt`. No torus
of the family can be larger than this regular region, which makes it a cheap way to see in
advance which members and mode mixes are worth continuing.

`breakdown_analysis.jl` combines the scans, the computed families and new orbit sweeps into the
"limits" figures: extent, frequency drift, and Fourier modes needed vs size.

## 10. `bin/param` reference

```
bin/param <input.csv> [options]
```

Continuation:

| option | default | meaning |
|---|---|---|
| `--domega d1 d2` | 0 0 | continuation direction in ω per unit ε. **Without it, ε does not change the problem** and every step recomputes the same torus |
| `--eps-max E` | 0.01 | stop beyond E (use a negative value for inward runs) |
| `--max-steps N` | 10 | maximum continuation steps |
| `--tol-floor T` | 1e-10 | accept a stalled Newton iteration whose best sup-norm error is below T LU |

Method switches, for experiments. The defaults are what produced the results:

| option | effect |
|---|---|
| `--free-omega` | Newton also corrects ω, and fixes the average normal correction instead (the energy is fixed by the section, so ω is otherwise an input) |
| `--no-filter` | do not zero the upper half of the spectrum after each Newton step (low-pass filtering is on by default) |
| `--no-real` | keep the imaginary part of K |
| `--fixed-frame` | do not rebuild the torus-adapted local coordinates after each accepted torus |
| `--tol-coho t` | relative threshold below which Fourier coefficients of the cohomological solutions are set to zero (default: the grid tolerance) |

Tools. Each one is the second argument and replaces the continuation:

| mode | what it does |
|---|---|
| `--test-twist [Case]` | Newton on the integrable twist map (§2) |
| `--fdcheck` | analytic DF vs finite differences at the first grid point |
| `--orbit niter q1 q2 p1 p2` | iterate the section map from z. Prints each iterate; the first line (`# Dfz …`) is DF(z), row-major |
| `--eval` | print F(K) and DF(K) at every grid point of the input torus: `l F0..F3 DF00..DF33` |

The input file of the tool modes supplies H, μ and the grid. `--orbit` and `--eval` are how the
Julia tools use the C++ map, so Julia and C++ always share the same integrator and section.

## 11. Troubleshooting

**The run stops after ten tori, or at ε = 0.01.** You left `--max-steps` or `--eps-max` at their
defaults.

**Every torus has the same ω.** `--domega` is missing.

**The first step fails (`Newton method does not converge` at ε = 0).** The start torus is not
close enough to invariant. Check the first `Error of invariance` line: an orbit-fitted start
should be ~1e-11. Larger values mean a poor fit (an irregular orbit, or too few Fourier modes).
Rerun `family_pipeline.jl` with a smaller amplitude or another mode mix.

**Steps shrink endlessly with errors just under `--tol-floor`.** The torus is truncated on its
grid, so each step barely passes. The proactive refinement (§5) normally catches this. If it
does not, restart from the last torus on a finer grid, or relax `--tol-floor` to 1e-10.

**Errors stall at ~1e-12 and never reach `tolinva` = 1e-12.** This is the noise floor of the map
evaluation, and it is expected (§5). Do not lower `--tol-floor` below the floor, or no torus will
ever be accepted.

**The grid jumps to 128×128 early on an inward branch.** This is a known effect of the older
refinement rule (see `results/families/nrho_id53/README.md`, `inward_1`). The current code only
refines on a stall within 100× of the floor. Restart from a resampled 32×32 torus if needed.

**Restarted run behaves differently from the original.** Check field 1 of the new start file
(§7.3). Also remember that ε restarted at 0.

**macOS: `clang: error: unsupported option '-fopenmp'`.** Build with Homebrew GCC (§2).

**`map failures` > 0.** An orbit failed to return to the section (escape or impact). Near the
edge of the regular region this is expected. On the start torus it means the torus is too large
for that member.

## 12. Map of the source code

`src/cpp/param.cc` (≈2200 lines) holds the solver:

| where | what |
|---|---|
| `main` (line ~308) | option parsing; reading the input; tool modes (`--orbit`, `--test-twist`, `--eval`, `--fdcheck`); the continuation loop (~670–950): predictor, Newton loop with stall detection and anisotropic refinement, step control, writing `output_torus`, rebuilding the local frame |
| `build_local_frame` (~168) | torus-adapted coordinates z = z_cen + M ζ, built from the first harmonics of the current torus, so that the torus looks like a product of round circles of size ~1. All tolerances are converted with `zscale`, and all printed errors are in physical LU |
| `kam_torus` (~1022) | one quasi-Newton step of the parameterization method: invariance error, tangent frame L = DK, normal frame N, torsion T, the two cohomological equations solved in Fourier space, update of K, tail check, low-pass filter |
| `map_CR3BP` (~1729) | the Poincaré map F and its differential DF, by integrating the flow and variational equations with the RK78 integrator of the Barcelona library up to the next q₃ = 0 upward crossing |
| `sform_CR3BP`, `gform_CR3BP`, `normal0_CR3BP` | symplectic form and metric in local coordinates, and the initial normal frame |
| `map_twist` (~1989) | integrable twist map for `--test-twist` |

`src/cpp/` also holds the Barcelona RTBP library: `rtbphp` (CR3BP Hamiltonian vector field and
variational equations), `seccp` (Poincaré sections), `rk78vp` (Runge–Kutta–Fehlberg 7(8): 13 stages, nodes 2/27, 1/9, …, weights 41/840, …),
`fluxvp`, and the grid/FFT and matrix classes in `headers/`. `docs/PROVENANCE.md` records who
wrote what.

`src/julia/`: see [`src/julia/README.md`](../src/julia/README.md).

## 13. Building the paper

```bash
cd docs/paper && make       # pdflatex + bibtex + 2 passes → paper.pdf
```

`make quick` runs a single pass without the bibliography. Figures are read from
`docs/paper/figures/`.
