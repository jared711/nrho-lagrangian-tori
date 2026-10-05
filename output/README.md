# output/

Not used by the current workflow. `bin/param` writes its converged tori, `output_torus<epsilon>`
(six decimals), into the directory it is run from; run each continuation branch in its own
directory (see `docs/STATUS_<date>.md`). Results worth keeping go to `results/families/`.

File format: 12 header lines (toltail, tolinva, tolinte, ω₁, ω₂, ε, H, μ, n₁, n₂, Δε, 1), then one
line per grid point: `i1 i2 q1 q2 p1 p2` (0-based grid indices, physical section coordinates).
This directory is gitignored except for this file.
