# NRHO Id 37 — tori with mode mix (0.3, 1) (outward, short)

Second slice of the Id 37 torus family (the main one, `../inward` and `../outward`, uses the mix
(1, 1)). Start torus from `family_pipeline.jl 37 <dir> 0.3 1`, continued outward (commit `302fbdd`,
tol-floor 1e-10). Only 22 tori, radius up to 7.31 km: the branch stopped early (steps rejected
on a 32×32 grid at the tolerance), well short of the ~16.7 km regular extent found by
`stability_edge.jl` in this direction.
