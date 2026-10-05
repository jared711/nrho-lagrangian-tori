# Outward extensions (2026-10-05)

Attempts to continue the outward branch of this family to larger tori, each restarted from the last
torus of the previous piece (ε relative to that torus; `src/julia/restart_from.jl`):

* `outward_ext1`: anisotropic grid refinement (code `b8d5e59`), tol-floor 1e-11.
* `outward_ext2`: relaxed tol-floor 1e-10 (~2.4 cm) and grid cap 262144 points (`f2b677c`).
* `outward_ext3`: as ext2, refinement only when the best torus is truncated (`302fbdd`).

Largest radius on the section: 12.28 → 12.4 km. The gain is small because the family has reached the edge
of the regular (KAM) region around this NRHO (`results/families/stability_edge_scan.txt`): beyond
it orbits become chaotic or escape, and the tori need rapidly growing numbers of Fourier modes.
To limit repository size, every piece keeps its full run log (all accepted tori are listed there with
their errors) but only representative tori: only the largest torus of ext3. The start files are
gzip-compressed restart inputs.
