# Julia tools

Run from the repository root with `julia --project=.`.

## Used for the current results

| file | purpose |
|---|---|
| `family_pipeline.jl` | `family_pipeline.jl <Id> <outdir>`: fixed point of the section map for NRHO `<Id>`, linear rotation numbers, a regular orbit-fitted 32×32 start torus, continuation direction Δω = ρ − ρ_lin; writes `start.csv`, `setup.txt` (commands) and `results/families/nrho_id<Id>/nrho_fixed_point.txt` |
| `pmap_tools.jl` | library: drives `bin/param --orbit/--eval` (fixed points, eigenframe, NAFF rotation numbers, orbit fitting, frequency map, map smoothness test, dense collocation Newton used as a reference solver) |
| `family_analysis.jl` | `family_analysis.jl <prefix> inward=<dir[,dir]> outward=<dir[,dir]>`: family table (CSV) and figures (frequencies, size, error, section views); env `NRHO_ID`, `RHO_LINEAR` |
| `qpo_visualization.jl` | `qpo_visualization.jl <prefix> <nrevs> <torus files>`: integrates QPOs from tori; 3D view, distance to the NRHO, apolune cross-sections; env `NRHO_ID` |
| `dynamics.jl`, `util.jl` | CR3BP vector field and STM in the Barcelona-convention Hamiltonian form; coordinate conversions (`rv2pq`, `pq2rv`) |
| `parameterization-method.jl` | `P()` (Poincaré map of the 6D flow with its differential) and section coordinate maps `γ`, `γ⁻¹`; the rest is an August 2026 partial Julia port of `kam_torus` (step 1 only), not used |

## Older

| file | purpose |
|---|---|
| `main.jl` | 2024 script: integrates NRHO Id 72, builds a linear (eigenvector) torus with the linear rotation numbers and writes `data/initial_conditions/approxQPO.csv`; fixed 2026-10-05, superseded by `family_pipeline.jl` |
| `test_kam_torus.jl` | August 2026 test of the partial Julia port |
