# Inward continuation toward the NRHO (continues `../inward`)

Start: the last torus of `../inward` (ε = −0.896 there), converted to an input file. ε here is
relative to that torus: ω(ε) = ω_start + ε·Δω with the same Δω = (−2.3333e-4, −4.6605e-4).
Cumulative ε = −0.896 + ε. The run stopped at ε = −0.0837 (cumulative −0.980, i.e. the
amplitude-frequency detuning down to ~2% of the original start torus) when the continuation
step fell below its minimum on the 32×32 grid (sup errors approaching 1e-11; this run predates
the automatic grid refinement of `ab4d59e`). Code: `9c9c77f`.

    bin/param start.csv --domega -2.3333e-4 -4.6605e-4 --eps-max -0.094 --max-steps 400 --tol-floor 1e-11
