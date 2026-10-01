# Family of Lagrangian tori around NRHO Id 53

Saturn–Enceladus CR3BP (μ = 1.901109735892602e-7), section q3 = 0 (upward crossing), member Id 53 of
`data/initial_conditions/NRHO_L2.csv`, energy H = −1.50001738952906 (Barcelona convention).
Perilune 17.1 km above the 252.1 km mean radius of Enceladus. Linear rotation numbers of the fixed
point (`nrho_fixed_point.txt`): ρ_lin = (0.33909577956895687, 0.1521834704601715).

Start torus: orbit-fitted (`pmap_tools.jl`) torus of amplitudes (4e-5, 4e-5) LU,
ρ = (0.3393830664659224, 0.15261138018692144), 32×32, error 8.1e-12 → 6.9e-13 LU after one Newton step.
Continuation in the rotation vector along Δω = ρ − ρ_lin = (2.8729e-4, 4.2791e-4) per unit ε
(ε = −1 is the NRHO):

    bin/param start.csv --domega 2.8729e-4 4.2791e-4 --eps-max <limit> --max-steps 3000 --tol-floor 1e-11

Each piece restarts from the last torus of the previous piece; ε inside a piece is relative to its
start torus. Cumulative ε of the pieces:

| piece | cumulative ε | notes |
|---|---|---|
| inward_1 | 0 → −0.695 | code `ab4d59e` (refinement on stall without the 100× guard: refined to 128×128 after predictor failures) |
| inward_2 | −0.695 → −0.860 | restarted on 32×32 (resampled), code `ab0f…`/guarded refinement (`<100×tol`), refined again to 128×128 |
| inward_3 | −0.860 → −0.932 | OpenMP build `9eebd24`; stopped when the step fell below its minimum (errors ~1e-11) |
| outward_1 | 0 → 2.650 | code `ab4d59e` |
| outward_2 | 2.650 → 3.271 | guarded refinement, 64×64 → 128×128 |
| outward_3 | 3.271 → 3.486 | OpenMP build `9eebd24`; stopped at errors ~1e-11 on 128×128 |

The accepted tori satisfy the sup-norm error of invariance < 1e-11 LU (~2.4 mm) by construction.
