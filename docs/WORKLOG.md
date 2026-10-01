# Work log

Running log of the solver work, newest last. Commit hashes refer to this repo. See
`PROVENANCE.md` for who wrote what.

## 2026-10-01 (Claude Opus 5.5 session with Jared)

### Bugs found and fixed
| Commit | Fix |
|---|---|
| `80585d6` | STM read transposed in `map_CR3BP` (column-major storage via `vr1()`). |
| `c45f46c` | `get_dp3()` returned ∇(p3²) instead of ∇p3, and read `z[4]` out of bounds; `gform_CR3BP()` had a typo. |
| `5a16052` | Leftover angle lift θ+ω in K(θ+ω); the invariance error was O(1) wrong (2.59). |
| `6ee5174` | `seccp()` failures were ignored. |
| `d9ab1ca` | Continuation loop used Case 1 (constant N₀), whose frame is singular for tori around an elliptic point. |
| not yet | `main.jl` writes ω = ρ/T₀ (rad/TU); `param` needs rotation numbers ρ/2π. |
| not yet | ε does not enter `map_CR3BP`, so the continuation loop re-solves the same torus. |

### Verified
* Map Jacobian vs finite differences: 5.7e-7 relative (`--fdcheck`).
* Map symplectic for the standard Ω to 5e-10.
* Fixed point (NRHO Id 72) by Newton: quadratic (3.9e-7, 2.5e-11, 2.7e-14). It is elliptic-elliptic, with ρ = (0.343616643, 0.173806970).
* `kam_torus` converges quadratically on an integrable twist map (`--test-twist`, `b32ab44`): 2.2e-4, 1.5e-6, 8.0e-11, 5.7e-13.

### Open problem: CR3BP Newton diverges from a very good guess
* Orbit-fitted guess (`pmap_tools.jl`, amplitude 1.2e-5 LU): invariance error 3.8e-12. One Newton step raises it to ~1e-8.
* ⟨T⟩ is nearly singular (det 3.6e-7 in local coordinates) and non-symmetric.
* Torus frequencies are near a 1:2 resonance: ρ₁ − 2ρ₂ = −0.0031.
* Hypotheses being tested: (a) integration noise of the map (~1e-12 absolute = ~1e-7 local) is comparable to the weak twist; (b) this orbit/amplitude is too close to resonance. Next: measure the noise, try larger amplitudes and tighter tolerances, scan the NRHO family for better-conditioned elliptic orbits.

### Update (later on 2026-10-01)
* Frequency map (`bc66033`): 1:2 resonance zone for a₁ ≳ 2e-5 LU; chaos for a₁ ≳ 8e-5; a regular region at a₁ ≲ 5e-6, a₂ ∈ [2e-5, 8e-5]. The first torus tried (1e-5, 1e-5) sits at the edge of the resonance zone.
* Map smoothness (`bc66033`): the map is smooth to ≲1e-13, so there is no noise floor. |D²P| ≈ 2.2e5, i.e. a nonlinearity scale of ~1 km.
* Tried; does not fix the quasi-Newton divergence: local coordinates (`14f7e02`), identity metric (`291b589`, which does improve the frame check 100×), free-frequency step (`5587d1c`, validated exactly on the twist map).
* Dense collocation Newton (`115274c`) **converges** from the same guess: 4.8e-10 → 4.9e-12, plateau 5e-12 (32×32). So the setup is consistent, and the quasi-Newton basin is the problem.
* Hybrid test: dense-Newton-cleaned torus (E = 2.3e-11) → kam_torus still diverges, both fixed-ω (→ 1.1e-8) and free-ω (→ 4.5e-9). Linearized residual 100–300× E.
* The torus is Lagrangian to its error level (W/(2πs)² ~ 4e-7, mean ~1e-13), so a non-Lagrangian guess is ruled out.
* Next: high-level review by three agents (literature on flow-map vs section-map methods; fresh code audit of kam_torus; numerical strategy).

### Update: quasi-Newton fixed; first tori and a short family (2026-10-01, evening)
* Correction: the stall was **not** a basin problem (as `115274c` and the entry above claimed). Agents found three bugs: one-sided Nyquist modes (`c76550b`), Newton step before the restart torus was saved (`cf9f7be`), and an under-resolved normal frame on eccentric tori (fixed by first-harmonic coordinates, `c4a0a83`, and pointwise N(θ+ω), `28700b7`). The final culprit was the absolute small-coefficient threshold in `cohomological()`, which zeroed ~85% of η_N for a torus of radius ~1e-5 (fixed by unit-circle coordinates `2f527c3` and a relative threshold `73ea3f0`).
* Map noise floor: section tolerance 1e-14 (`04dcb00`). Floor ≈ 1.2e-11 LU (~3 mm) on a ~10 km torus. It is not set by the RK78 tolerance, `-ffast-math`, or grid size (it gets worse on finer grids: small divisors amplify the noise).
* Stall detection (`491408d`) and Haro–Mondelo low-pass filter (`03032b6`): Newton converges 2.6e-10 → 1.2e-11 and stays there. **First accepted CR3BP torus** (NRHO Id 97).
* Continuation in ω (`912714b`): 5 tori accepted up to ε = 0.125 (|Δω| ~ 6e-5), then stuck as the floor rises toward tol-floor. Refining the grid on stall made it worse (not committed).
* Open: the source of the ~1e-11 floor; continuation step size; nearest resonance along the path is (−1, 7) at 3e-3.
