# Tests

The checks used for this project are built into `bin/param`:

* `bin/param <input.csv> --test-twist 2` — integrable 4D twist map with known tori: Newton must
  converge quadratically (2.1e-4, 1.3e-6, 7.1e-11, 1.1e-16 with the default input grid).
* `bin/param <input.csv> --fdcheck` — Jacobian of the CR3BP section map vs central finite
  differences (relative agreement ~6e-7, limited by the finite-difference step).
* `src/julia/pmap_tools.jl`: `fixed_point` (quadratic Newton on the NRHO), `map_noise`
  (smoothness of the map), `collocation_newton` (dense reference Newton step).

`test.cpp` is a leftover from the first commit (2024) and is not part of the build.
