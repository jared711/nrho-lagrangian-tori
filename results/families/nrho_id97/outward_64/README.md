# Outward continuation with automatic grid refinement (continues `../outward_32`)

Start: the last torus of `../outward_32` (ε = 0.570 there), on 32×32; the first Newton loop
stalls at 9.9e-12 and refines the grid to 64×64. ε here is relative to that torus (cumulative
ε = 0.570 + ε); 14 tori up to ε = 2.03 (cumulative 2.60), sup errors 3.9e-13 .. 6.6e-13 LU.
Code: `ab4d59e`.

    bin/param start.csv --domega -2.3333e-4 -4.6605e-4 --eps-max 2.43 --max-steps 1000 --tol-floor 1e-11
