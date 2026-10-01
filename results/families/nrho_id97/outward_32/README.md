# Outward continuation on a 32×32 grid

Start: the same start torus as `../inward` (column 11, the step, set to +0.005). 16 tori up to
ε = 0.570; then the sup error reaches ~1e-11 (truncation on 32×32) and the continuation stops.
Continued on 64×64 in `../outward_64`. Code: `9c9c77f`.

    bin/param start.csv --domega -2.3333e-4 -4.6605e-4 --eps-max 3 --max-steps 1000 --tol-floor 1e-11
