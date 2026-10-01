# Family of Lagrangian tori around NRHO Id 23

Set up with `src/julia/family_pipeline.jl 23` (commit `b9fa65e`); continued with `bin/param` at commit `8e79f0e` (OpenMP, -O2, refinement on stall within 100x of tol-floor), 4 threads per branch. Perilune altitude km above the 252.1 km mean radius of Enceladus.

`setup.txt`: rotation numbers, continuation direction and the commands used. `inward/` and `outward/`: start torus, converged tori and run log (gzip-compressed). Both branches ended when the continuation step fell below its minimum.
