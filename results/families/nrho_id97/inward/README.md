# Family of Lagrangian tori around NRHO Id 97 (inward branch)

Saturn–Enceladus CR3BP (μ = 1.901109735892602e-7), section q3 = 0 (upward crossing), energy
H = −1.50001746408344 (Barcelona convention), member Id 97 of `data/initial_conditions/NRHO_L2.csv`.

* `start.csv`: input torus. It is an orbit-fitted torus (Julia `pmap_tools.jl`: `fixed_point`,
  `eigen_frame`, `refined_rotation_numbers`, `fit_torus`) of amplitudes (4e-5, 4e-5) LU along the
  two elliptic eigendirections. ρ = (0.35372425999973445, 0.19383759358522312), 32×32 grid. The
  continuation step field (column 11) is −0.01.
* `output_torus<ε>`: converged tori (format written by `param`: 12 header lines with toltail,
  tolinva, tolinte, ω₁, ω₂, ε, H, μ, n₁, n₂, Δε, 1, then one line per grid point:
  `i1 i2 q1 q2 p1 p2`, physical coordinates, 0-based indices).
* `run.log`: full output of the run, including each Newton history.

Command (binary built from commit `9c9c77f`):

    bin/param start.csv --domega -2.3333e-4 -4.6605e-4 --eps-max -0.9 --max-steps 400 --tol-floor 1e-11

ω(ε) = ρ + ε·Δω, with Δω = ρ − ρ_linear (the amplitude-frequency detuning of the start torus), so
ε = −1 is the NRHO itself. 43 tori from ε = 0 to −0.896; sup error of invariance: median 4.0e-13 LU,
worst 7.2e-12 LU.
