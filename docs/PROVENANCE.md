# Code provenance

Who wrote what in this repository. It keeps the credit straight for the paper,
and it shows which code was validated by its original authors and which code is
new. Keep it up to date: every commit that changes code should be traceable to
a row here.

Sources:
* **BCN-KAM**: Àlex Haro (UB) and Alejandro Luque (UB). The parameterization-method
  KAM torus code (header of `param.cc`: "Last version October 21, 2014").
* **BCN-RTBP**: Àlex Haro and Josep-Maria Mondelo (UAB). The CR3BP vector field,
  integrator and Poincaré-section library (`rtbphp.c`: "Copyright (C) 2016 Àlex
  Haro, Josep-Maria Mondelo", LGPL v3 or later; Catalan comments throughout).
* **JB-2024**: Jared Blanchard, May–July 2024 (commits `132a7ea`..`ca0df7c`).
* **AI-2026-08**: Claude Sonnet 4.5 sessions, 2026-08-04/05 (commits `f160d04`..`6d53f3c`,
  co-authored trailers).
* **AI-2026-10**: Claude Opus 5.5 session, 2026-10-01 onward (co-authored trailers).

The first commit (`132a7ea`) already mixes BCN-KAM and JB-2024 code in `param.cc`, so
the split inside that file comes from the code itself (comments, function names), not
from git history.

## C++ (`src/cpp/`)

| File | Origin | Later changes |
|---|---|---|
| `headers/grid.h`, `headers/matrix.h`, `headers/complex.h` | BCN-KAM | JB-2024: debug prints added to the `matrix` copy constructor and `fft_F` (`ca0df7c`). AI-2026-10: Nyquist modes zeroed in `deriva`, `shift`, `cohomological`; `cohomological` threshold made relative (`grid.h`) |
| `rtbphp.c`, `seccp.c`, `fluxvp.c`, `rk78vp.c`, `campvp.c`, `scread.c`, `vbprintf.c`, `headers/{rtbphp,seccp,fluxvp,rk78vp,campvp,scread,vbprintf,utils}.h` | BCN-RTBP | JB-2024: English translations of the Catalan comments, argument documentation in `seccp.c`, reformatting. The "RTBP" integrator settings in `fluxvp.c` (`fluxvp_tol=4e-14`, etc.) were already in the first commit, so they are BCN-RTBP |
| `param.cc`: `kam_torus()`, `realloc_torus()`, the continuation/Newton driver in `main()`, and the standard and Froeschle maps | BCN-KAM | see below |
| `param.cc`: `nu()`, `get_p3()`, `get_dp3()`, `map_CR3BP()`, `sform_CR3BP()`, `gform_CR3BP()`, `normal0_CR3BP()`, `wrtf()`, `state2ham()`, the `DTOR/DMAP/NPAR` settings for the CR3BP, and the input-file reading | JB-2024 ("functions created by Jared Blanchard May, 2024") | see below |
| `Makefile` | JB-2024 | AI-2026-08: new directory layout (`1668c0c`) |

Changes to BCN-KAM code inside `param.cc`:
* JB-2024: commented out the angle lift `z[i] += index[i]/nn[i]` where K is evaluated (dated 6/27/24), and removed the identity from `LR` ("Alex 5/3 get rid of the val1", i.e. advice from Àlex Haro).
* AI-2026-08 (`c52cac6`): commented out the second `kam_torus()` call and added `MAX_CONT_STEPS`/`MAX_EPSILON`. **The diagnosis behind this change was wrong.** Repeating a Newton step is correct; the error grew because of the bugs fixed in AI-2026-10.
* AI-2026-10 (`5a16052`): removed the remaining angle lift on K(θ+ω).
* AI-2026-10: removed the unconditional `kam_torus()` call before the continuation loop (found by a code-audit agent).
* AI-2026-10: `realloc_torus()` drops Nyquist coefficients when refining the grid (found by a code-audit agent).
* AI-2026-10 (`6ee5174`): `kam_torus()` returns −1 when the Poincaré map failed.
* AI-2026-10 (`f83e8f0`): `kam_torus()` prints diagnostics after each step: the residual of the linearized equation, the Lagrangian defect, and the symplectic-frame check. `3c7ad63`: also the averaged torsion ⟨T⟩.
* AI-2026-10: `map_twist()`, `gform_identity()` and the `--test-twist` mode, an integrable twist map that validates `kam_torus()` on Cartesian (non-lifted) tori.
* AI-2026-10: the Newton loop in `main()` uses Case 2 (metric frame) in place of Case 1.
* AI-2026-10: both `kam_torus()` calls for the CR3BP use `gform_identity()` (the Euclidean metric in local coordinates) instead of `gform_CR3BP()`.
* AI-2026-10: `--free-omega` option in `kam_torus()`: the Newton step corrects ω and fixes the average normal correction, instead of fixing ω and inverting ⟨T⟩.
* AI-2026-10: `--eval` mode in `main()` (F and DF at every grid point, for the Julia collocation solver).

Changes to JB-2024 code inside `param.cc` (all AI-2026-10):
* `80585d6`: STM index order in `map_CR3BP()`.
* `c45f46c`: `get_dp3()` chain rule and index; `gform_CR3BP()` now reuses it.
* `6ee5174`: check of the `seccp()` return value.
* AI-2026-10: section-crossing tolerance `tolJM` in `map_CR3BP()` 1e-12 → 1e-14.
* `4b6c47a`: `--fdcheck` and `--orbit` modes in `main()`.
* AI-2026-10: local coordinates ζ = (z − z_c)/s for the CR3BP torus (`to_physical()`, `zcen`, `zscale`); `map_CR3BP()` and `gform_CR3BP()` take ζ.
* AI-2026-10: local coordinates generalized to z = z_c + M ζ, with M from the symplectically normalized first harmonics of the input torus, rescaled so the torus circles have radius ~1 (`Mloc`, `Minv`, `Omega_loc`, `invert4()`); `sform_CR3BP()` returns Ω_loc = MᵀΩM (suggested by the code-audit and literature agents).
* AI-2026-10: in `kam_torus()` (Case 2), N(θ+ω) is computed pointwise from L(θ+ω) with the frame formula instead of Fourier-shifting N (suggested by the code-audit agent).

## Julia (`src/julia/`)

| File | Origin | Later changes |
|---|---|---|
| `main.jl`, `dynamics.jl`, `parameterization-method.jl` (first 66 lines), `util.jl` | JB-2024 (uses Jared's `ThreeBodyProblem.jl` package) | AI-2026-08: `Vern9` in place of `TsitPap8`, `global uidx`, the `terminate!` import (`f160d04`), new file paths (`7e3f235`) |
| `parameterization-method.jl` (KAM STEP 1 port), `test_kam_torus.jl`, `README.md` | AI-2026-08 (`6d53f3c`) | |
| `pmap_tools.jl` (fixed points, NAFF frequency analysis, torus fitting, frequency map, dense collocation Newton) | AI-2026-10 | |

## Data (`data/`)

* `initial_conditions/NRHO_L1.csv`, `NRHO_L2.csv`: JB-2024, added in `d202b3d`. These are Saturn–Enceladus halo/NRHO initial conditions (μ = 1.901109735892602e-7). TODO (Jared): record the original source.
* `initial_conditions/approxQPO*.csv`, `halo.csv`, `T₀.csv`: generated by JB-2024 `main.jl`. The rotation numbers in `approxQPO.csv` are in rad/TU (ρ/T₀); `param` needs ρ/2π.
* `config/*.txt`: from the first commit; origin unclear (probably BCN-KAM example inputs).

## Documentation

* `docs/paper/`, `docs/analysis/`, `docs/STATUS.md`, `docs/JULIA_IMPLEMENTATION_STATUS.md`, and the untracked `docs/REORGANIZATION_COMPLETE.md` and `docs/VERIFICATION_TESTS.md`: AI-2026-08. `docs/paper/references.bib` contains fabricated and incorrect entries (checked 2026-10-01) and must be rebuilt before use.
