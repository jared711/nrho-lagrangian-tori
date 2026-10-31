# Companion notes to the paper

These notes explain the harder ideas behind `paper.tex` at more length than a paper can: where the
equations come from, why each ingredient is there, and how the pieces fit together. They follow
the order of the paper. Notation matches the paper. For the code, see
[`docs/TUTORIAL.md`](../TUTORIAL.md).

Contents

1. [The chain of reductions](#1-the-chain-of-reductions)
2. [The CR3BP in Hamiltonian form](#2-the-cr3bp-in-hamiltonian-form)
3. [The Poincaré map at fixed energy](#3-the-poincaré-map-at-fixed-energy)
4. [NRHOs as elliptic fixed points](#4-nrhos-as-elliptic-fixed-points)
5. [Invariant tori, and why they are Lagrangian](#5-invariant-tori-and-why-they-are-lagrangian)
6. [KAM theory in one picture](#6-kam-theory-in-one-picture)
7. [The parameterization method step by step](#7-the-parameterization-method-step-by-step)
8. [What the a-posteriori theorem buys](#8-what-the-a-posteriori-theorem-buys)
9. [Why each practical ingredient is needed](#9-why-each-practical-ingredient-is-needed)
10. [Continuation in the rotation vector](#10-continuation-in-the-rotation-vector)
11. [Initial guesses: NAFF, frequency drift and orbit fitting](#11-initial-guesses-naff-frequency-drift-and-orbit-fitting)
12. [Breakdown of tori and the limits of the families](#12-breakdown-of-tori-and-the-limits-of-the-families)
13. [How the other methods in the literature work](#13-how-the-other-methods-in-the-literature-work)
14. [Glossary](#14-glossary)

---

## 1. The chain of reductions

The whole paper rests on a sequence of reductions, each removing one dimension or one symmetry:

| object | dimension | what removed it |
|---|---|---|
| CR3BP flow | 6 (q, p ∈ R³) | — |
| energy level H = h | 5 | H is conserved |
| Poincaré section q₃ = 0 in that level | 4 | one crossing per revolution replaces time |
| NRHO | periodic orbit → **fixed point** z\* of F | |
| QPO family around it | 3-tori of the flow → **2-tori** of F | |
| torus | parameterization K: T² → R⁴ | the unknown of the invariance equation |

Each reduction keeps the Hamiltonian structure: F is a symplectic map of R⁴. Everything
afterward, from Lagrangian tori and the symplectic frame to the torsion and the cohomological
equations with constant coefficients, depends on that structure.

## 2. The CR3BP in Hamiltonian form

In the rotating frame with positions q and momenta p,

$$H(q,p) = \tfrac12|p|^2 + q_2p_1 - q_1p_2 - \frac{1-\mu}{r_1} - \frac{\mu}{r_2}.$$

The terms $q_2p_1 - q_1p_2$ come from the rotation of the frame, which is why p is not the
velocity: $\dot q_1 = p_1 + q_2$, $\dot q_2 = p_2 - q_1$, $\dot q_3 = p_3$. The code uses the
Barcelona convention, in which the primaries sit at $q_1 = \mu$ (Saturn) and $q_1 = \mu-1$
(Enceladus). This is the standard rotating frame rotated by 180° about z. The two are related by
$H = -C/2$, where C is the Jacobi constant.

The Hamiltonian form matters because it makes the flow **symplectic**. It preserves the 2-form
$\Omega = \sum_i dq_i\wedge dp_i$, which measures oriented areas in each (qᵢ, pᵢ) plane. In
matrix form, a map is symplectic when its Jacobian M satisfies $M^\top\Omega M = \Omega$.
Symplecticity is a much stronger condition than preserving volume. It is the property that makes
KAM theory apply and that the parameterization method exploits.

## 3. The Poincaré map at fixed energy

**Dimension count.** Fix $H = h$ (5-D level) and cut it with $q_3 = 0$, keeping upward crossings
($p_3 > 0$). A point of the section is determined by $z = (q_1,q_2,p_1,p_2)$, and $p_3$ follows
from $H = h$ (eq. (3) of the paper, where the positive root is the upward crossing). F maps z to
the next upward crossing.

**Why F is symplectic in z.** The first-return map of a Hamiltonian flow, restricted to an energy
level and a transversal section, preserves the restriction of Ω to the section. On
$\{q_3 = 0, H = h\}$ that restriction is $dq_1\wedge dp_1 + dq_2\wedge dp_2$, because $dq_3 = 0$
there. So F is symplectic for the *standard* form in z, even though p₃ varies along the section.
This is why the code can use the constant matrix Ω.

**The differential.** Integrating the variational equations gives the state transition matrix Φ
over the return time τ(z). But τ depends on z: a perturbed orbit reaches the section earlier or
later. The correction term in

$$DP = \Phi + f(x_\tau)\,D\tau^\top,\qquad D\tau^\top = -\frac{\nabla g^\top\Phi}{\nabla g\cdot f(x_\tau)}$$

slides the perturbed endpoint along the flow, by $f\,\delta\tau$, until it lies back on the
section. $D\tau$ is chosen so that $\delta q_3 = 0$ at the endpoint. Then $DF = \pi\,DP\,D\nu$
adds the chain rule for p₃(z) at the start and drops the q₃, p₃ components at the end. The
`--fdcheck` test verifies exactly this formula against finite differences.

**Order-one entries.** The entries of DF are of order 27 (paper §3.3). The NRHO passes very
close to Enceladus, and the map strongly shears the section near perilune. This becomes
important in §9.

## 4. NRHOs as elliptic fixed points

A periodic orbit that crosses the section once per period is a fixed point $F(z^*) = z^*$. Since
F is symplectic, the eigenvalues of $DF(z^*)$ come in quadruples $\lambda, 1/\lambda,
\bar\lambda, 1/\bar\lambda$. In 4-D there are two pairs, and each can be

* **elliptic**: $\lambda = e^{\pm 2\pi i\rho}$ on the unit circle (rotation by angle $2\pi\rho$),
* **hyperbolic**: real, $\lambda$ and $1/\lambda$ with $|\lambda| > 1$ (a saddle), or
* in a complex quadruple off the circle (complex saddle).

The paper uses members that are **elliptic–elliptic**: both pairs lie on the circle, with linear
rotation numbers $\rho = (\rho_1,\rho_2)$. Near $z^*$, in the coordinates of the two eigenplanes
(the real and imaginary parts of the eigenvectors $V_1, V_2$), F is close to two independent
rotations by $\rho_1$ and $\rho_2$. Invariant tori are then approximately products of two circles,
one in each eigenplane, with radii given by the two mode amplitudes. That is the "mode mix"
$(c_1, c_2)$ of the paper: it chooses the ratio of the two radii.

Linear stability (elliptic) is necessary but not sufficient for nonlinear stability. KAM theory
(§6) supplies most of the rest. In 4-D maps, even KAM does not confine orbits forever (§6,
"no barriers").

## 5. Invariant tori, and why they are Lagrangian

A torus of F with rotation vector ω is the image of a map $K:\T^2\to\R^4$ that satisfies

$$F(K(\theta)) = K(\theta+\omega).$$

Iterating F on the torus is a rigid translation by ω in the angles θ. If ω is **non-resonant**
($k\cdot\omega\notin\mathbb Z$ for every integer $k\neq 0$), each orbit on the torus fills it
densely. Lifted back to the flow, the torus is a 3-torus with frequencies
$(1/T_{\text{return}}, \omega_1/T, \omega_2/T)$, roughly speaking, and the QPOs are orbits that
wind densely around it.

**Lagrangian** means that Ω vanishes on the tangent space:
$DK(\theta)^\top\,\Omega\,DK(\theta) = 0$. A 2-D surface in R⁴ carries zero symplectic area in
every tangent plane. Here is why every invariant torus with non-resonant ω has this property:

1. $DK^\top\Omega DK$ is an antisymmetric 2×2 matrix, so it is determined by one function $c(\theta)$.
2. Differentiate the invariance equation: $DF(K(\theta))\,DK(\theta) = DK(\theta+\omega)$. Since
   DF is symplectic, $c(\theta+\omega) = c(\theta)$. A function invariant under a dense
   translation is constant, so $c(\theta) = c$.
3. F is exact symplectic, which means $\Omega = d\alpha$ globally. The integral of $c$ over $\T^2$
   is the integral of the exact form $K^*\Omega = d(K^*\alpha)$ over a closed surface, so it is 0.
   Hence $c = 0$.

The Lagrangian property is what lets you build a frame (§7) whose tangent half is DK itself and
whose normal half is symplectically paired with it.

**Cartesian tori.** Around an elliptic fixed point the tori encircle the point. K is periodic in
both angles with no "lift" term like θ. The paper notes this because the textbook setting
(tori in action–angle coordinates) has $K(\theta) = (\theta + \ldots, I(\theta))$, which is not
periodic in θ.

## 6. KAM theory in one picture

**Integrable picture.** If F were integrable, there would be action–angle coordinates (I, φ) in
which $F(I,\varphi) = (I, \varphi + \omega(I))$. Every level set $I = \text{const}$ would be an
invariant torus. The frequency map $I\mapsto\omega(I)$ has a derivative, the **twist**
$\partial\omega/\partial I$. "Twist condition" means that this matrix is invertible, so that
nearby tori rotate at different rates.

**What survives a perturbation.** Kolmogorov–Arnold–Moser: for a small enough perturbation of a
twisting integrable map, the tori whose ω is **Diophantine** survive, slightly deformed.
Diophantine means

$$|k\cdot\omega - m| \ge \frac{\gamma}{|k|^{\tau}} \quad\text{for all } k\neq 0,\ m\in\mathbb Z.$$

This is a quantitative "very irrational". The surviving tori form a **Cantor set**: no interval
of tori survives, but the set has positive measure. Its relative measure tends to 1 as the
perturbation tends to 0. The gaps sit around resonant frequencies $k\cdot\omega\in\mathbb Z$.
Inside the gaps live resonant island chains (lower-dimensional tori, periodic orbits) and thin
chaotic layers.

**Why it applies near the NRHO.** Near an elliptic fixed point, the Birkhoff normal form makes F
approximately integrable. The distance from the fixed point plays the role of the perturbation
size, with the squared amplitudes as actions. The twist is the amplitude dependence of the
rotation numbers, which is why the paper needs ρ to change with torus size. So close to the NRHO
almost all tori survive (KAM). Farther out, the resonant gaps widen, the tori become less smooth,
and eventually they break up (§12).

**No barriers in 4-D.** In a 2-D map, an invariant circle separates the plane and traps orbits
inside it forever. In a 4-D map a 2-torus has codimension 2 and separates nothing, like a circle
in 3-D space. Orbits that start between tori can in principle wander away through the web of
resonances (Arnold diffusion), though usually extremely slowly. The paper's tori are invariant
objects, and orbits *on* them stay on them. Nearby off-torus orbits are long-term stable only in
a practical, not an absolute, sense.

**Two-parameter family.** In 4-D there are two actions, so the tori form a two-parameter
(Cantor) family indexed by $\omega\in\R^2$, or equivalently by the two mode amplitudes. A
continuation along a line in ω (§10) follows one path through this family.

## 7. The parameterization method step by step

Unknown: K (on an $n_1\times n_2$ grid), with ω fixed. Error:
$E(\theta) = F(K(\theta)) - K(\theta+\omega)$.

**Plain Newton** would solve $DF(K(\theta))\Delta K(\theta) - \Delta K(\theta+\omega) = -E(\theta)$
for ΔK. Discretized, this is a dense $4N\times4N$ linear system: the shift by ω couples every
grid point to every Fourier mode. That costs $O(N^3)$, and the system is also singular, because
$K(\theta + \phi)$ is a solution for any phase φ.

**The trick: change frame so the linear operator becomes (almost) constant.** At each θ build a
4×4 matrix $P(\theta) = (L(\theta)\ \ N(\theta))$:

* $L = DK$: two columns tangent to the torus.
* $N$: two normal columns, chosen so that $L^\top\Omega N = -I$ and $N^\top\Omega N = 0$ (up to
  O(E); the code prints $|L^\top\Omega N + I|$ as its frame check). This makes P symplectic. The paper's formula (6), with metric G, is one explicit way to
  build N from L.

In this frame the linearized map has a simple form. Write $DF(K(\theta))P(\theta) =
P(\theta+\omega)\Lambda(\theta)$. Then:

* **Tangent block = I.** Differentiating the invariance equation gives $DF\,L(\theta) =
  L(\theta+\omega)$ exactly (when E = 0). The map carries tangent vectors of the torus to tangent
  vectors, and the first two columns of Λ are $(I;0)$.
* **Lower-right block = I.** This follows from symplecticity and the Lagrangian property. The
  quantity $\Lambda^\top\Omega\Lambda = \Omega$ and the first block column being $(I;0)$ force the
  lower-right block to be $I$.
* **Upper-right block = T(θ), the torsion.** It is the one block not fixed by the geometry, and it
  measures how normal displacements are sheared along the torus. It is the twist of §6, seen
  pointwise.

So, up to O(E), $\Lambda = \begin{pmatrix} I & T\\ 0 & I\end{pmatrix}$. Write $\Delta K = P\xi$ with
$\xi = (\xi_L,\xi_N)$ and project the error with $P(\theta+\omega)^{-1}$, which is cheap because
P is symplectic: $P^\top\Omega P = \Omega$ gives $P^{-1} = \Omega^{-1}P^\top\Omega$, whose rows produce the paper's $\eta_L, \eta_N$.
The Newton equation becomes triangular:

$$\xi_N(\theta) - \xi_N(\theta+\omega) = \eta_N(\theta),$$
$$\xi_L(\theta) - \xi_L(\theta+\omega) = \eta_L(\theta) - T(\theta)\,\xi_N(\theta).$$

**Cohomological equations.** Each line has the form $u(\theta) - u(\theta+\omega) = v(\theta)$.
In Fourier space the shift is a multiplication, $\widehat{u(\cdot+\omega)}_k =
e^{2\pi i k\cdot\omega}\hat u_k$, so

$$\hat u_k = \frac{\hat v_k}{1 - e^{2\pi i k\cdot\omega}},\qquad k\neq 0.$$

The system is diagonal, so it costs one FFT each way. Two subtleties:

* **Small divisors.** $1 - e^{2\pi ik\cdot\omega}$ is small when $k\cdot\omega$ is close to an
  integer. The Diophantine condition bounds how small it can get, $\gtrsim \gamma/|k|^\tau$. Then
  u is a little less smooth than v, but still smooth. This is the "loss of derivatives" that the
  KAM proof controls. Numerically, it is why high-k noise is amplified (§9, filter).
* **Averages.** The k = 0 equation reads $0 = \hat v_0$. The equation is solvable only if v has
  zero average, and then u is determined up to a constant. For the normal equation the free
  constant is $\xi_N^0$, which is chosen to make the tangent equation solvable:
  $\langle \eta_L - T\xi_N\rangle = 0$, which gives $\xi_N^0 = \langle T\rangle^{-1}\langle\ldots\rangle$.
  This is where the **twist condition** $\det\langle T\rangle\neq 0$ enters. The free constant of
  the tangent equation, $\xi_L^0$, is a shift of the phase θ, which is the singular direction of
  plain Newton. Setting it to 0 fixes the phase.

**Why it is still quadratic.** The frame identities hold up to O(E). Dropping those terms adds an
error of $O(E^2)$ to the step, the same order as Newton's own neglected terms. The iteration
therefore converges quadratically: the twist-map test shows the error going roughly from e to e².

**Cost.** Every operation is pointwise on the grid (building P, multiplying by T) or diagonal in
Fourier space. So a step is $O(N\log N)$ plus N evaluations of F and DF, and in the CR3BP the map
evaluations dominate. Compare $O(N^3)$ and $O(N^2)$ memory for a dense Newton solve (§13).

**The frequency-adjusting variant** (`--free-omega`) sets $\xi_N^0 = 0$ and instead changes ω so
that the tangent equation is solvable. It needs no $\langle T\rangle^{-1}$, which is useful as a
test or when the twist is weak. The paper computes families at fixed ω.

## 8. What the a-posteriori theorem buys

"A-posteriori" means that the theorem takes an *approximate* torus as input, such as the one the
computer returned, and concludes that a true torus exists nearby. Roughly:

> If ω is Diophantine (γ, τ), and the approximate torus K has error ‖E‖ small compared with an
> explicit expression in γ, τ, ‖DK‖, ‖N‖, ‖⟨T⟩⁻¹‖ and the width of the analyticity strip, then
> there is a true invariant torus K\* with the same ω and ‖K\* − K‖ ≤ C‖E‖.

The proof runs the same quasi-Newton iteration as the code and shows that it converges. This is
why the method is called "KAM-inspired". Figueras–Haro–Luque evaluate the constants with interval
arithmetic to obtain computer-assisted proofs. The paper does not do this. But the theory tells
you which quantities to watch: small error, moderate ‖⟨T⟩⁻¹‖ (twist not degenerate), and Fourier
coefficients that decay (an analytic torus). These are exactly the diagnostics in the log.

## 9. Why each practical ingredient is needed

The tori are tiny ($10^{-5}$–$10^{-4}$ LU), sit at |z| ≈ 1, are very eccentric, and live under a
map with entries of order 27. Each ingredient addresses one consequence of that.

**Torus-adapted coordinates (the most important).** In physical coordinates every quantity mixes
the O(1) position of the torus with its O(10⁻⁵) size and its eccentric shape. Then:

* $L = DK$ is tiny and badly conditioned, and $L^\top GL$ in the frame formula has a condition
  number of order (eccentricity)²;
* absolute tolerances designed for O(1) tori are meaningless.

The fix is to change coordinates affinely, $z = z_c + M\zeta$. The columns of M are the
first-harmonic axes of the torus, normalized symplectically and rescaled, so that in ζ the torus
is close to two unit circles. This is the same idea as a symplectic Floquet normalization: work in
coordinates adapted to the linear dynamics, except that here they are taken from the torus
itself. The map becomes $\tilde F(\zeta) = M^{-1}(F(z_c+M\zeta) - z_c)$, everything is O(1), and
the Euclidean metric G = I is reasonable. Along the family the torus changes shape, so the
coordinates are rebuilt after every accepted torus.

**Scale-free thresholds.** A cut-off like "skip Fourier coefficients below 1e-15" means
different things for a torus of size 1 and one of size 10⁻⁵. Thresholds relative to the largest
coefficient behave the same way at every scale.

**Nyquist modes.** On an even grid, mode $k = n/2$ is its own alias, and multiplying it by
$e^{2\pi ik\omega}$ or by $ik$ breaks the conjugate symmetry that keeps K real. Zeroing it avoids
spurious imaginary parts.

**Shifted normal frame.** N contains $(L^\top GL)^{-1}$. That makes it less smooth than L, so its
Fourier series converges more slowly, and shifting N by Fourier interpolation loses accuracy.
Computing $N(\theta+\omega)$ from the shifted L is exact to the accuracy of L.

**Sup-norm error in physical units.** The ℓ¹ norm of all Fourier coefficients adds up the noise
of all N coefficients. It grows with the grid, so refining the grid made the error look worse.
The maximum over the grid in physical LU is grid-independent and means something physical
(1e-11 LU ≈ 2.4 mm).

**Low-pass filter.** DFT coefficients at high |k| have the largest relative error: aliasing, and
round-off relative to their small size. The small divisors then amplify exactly those modes.
Without the filter the iteration seems to converge and then diverges. Zeroing the upper half of
the spectrum after each step ($|k_j|\ge n_j/4$) removes the modes that cannot be trusted. The
price is that the grid must be twice as fine as the torus needs.

**Noise floor and stall detection.** F is computed by numerical integration with local
tolerance 4e-14, plus a section-crossing refinement. Its error is propagated through a map with
large derivatives, so F itself is only accurate to about 1e-13–1e-12 LU. Newton cannot push E
below the error of F. Beyond that point the step is driven by noise, and the small divisors
amplify it. So the code keeps the best torus, stops when the error stops falling, and accepts the
torus if it is below the floor $\epsilon_{\text{floor}}$ = 1e-11 LU.

**Refinement, and why anisotropic.** If Newton stalls above the floor and the error spectrum sits
just below the filter cutoff, the torus is truncated: it needs more modes. Doubling only the
angle that is truncated keeps the cost growing by 2× instead of 4×. Near the edge of the family
one angle needs far more modes than the other (127 vs 63 for the largest Id 51 torus), so this
saving makes a large difference there.

## 10. Continuation in the rotation vector

The energy is fixed by the section, so the free parameters of a torus are its two frequencies
(§6). The paper continues along a line

$$\omega(\varepsilon) = \omega_s + \varepsilon\,\Delta\omega,\qquad \Delta\omega = \rho_s - \rho_{\text{lin}},$$

where $\rho_s$ are the rotation numbers of the start torus and $\rho_{\text{lin}}$ those of the
NRHO. So ε = 0 is the start torus and ε = −1 is the NRHO. Tori shrink as ε → −1 and grow for
ε > 0. Since ω depends nearly linearly on the squared amplitudes near the NRHO, a line in ω is a
sensible path. Its direction is set by the start torus, and so by the mode mix (paper §6.5).

**But the family is a Cantor set: how can a line work?** A line in ω meets resonances densely.
In practice the method still converges because Newton only needs ω to be far enough from
*low-order* resonances at the grid's resolution. Resonances of order |k| above the retained
modes are invisible to the discretization, and their gaps are exponentially thin for analytic
tori. Near a low-order resonance the twist condition and the small divisors deteriorate,
continuation steps shrink, and the branch may stop. That is why the member scan reports the
distance to the nearest low-order resonance.

**Step control** (Haro–Mondelo). The predictor is a secant through the last two tori. If Newton
needs fewer than 4 iterations, the step grows; if more, it shrinks; if Newton fails, the step is
halved. The branch ends when the step falls below 1/1000 of its initial value.

## 11. Initial guesses: NAFF, frequency drift and orbit fitting

**Why orbits make good start tori.** An orbit started on a KAM torus stays on it and fills it
densely. Its iterates $z_m = F^m(z_0)$ sit at angles $\theta = m\rho$ on the torus (with a
suitable origin), so fitting $z_m \approx K(m\rho)$ by least squares recovers K.

**NAFF (Laskar).** To get ρ from an orbit, project the iterates onto the complex eigen-coordinate
of each elliptic mode and find the dominant frequency of that complex signal. A linear fit of the
unwrapped angle gives a first guess. Maximizing the Hann-windowed Fourier amplitude near that
guess refines it to much better than 1/(number of iterates), because the window suppresses
leakage.

**Frequency drift as a regularity test.** Compute ρ from the first and second halves of the
orbit separately. On a KAM torus both halves have the same frequency, and the drift is limited
only by the accuracy of NAFF (1e-11–1e-9 here). On a chaotic orbit, or one trapped in a resonance
layer, the frequency wanders, and the drift jumps by orders of magnitude. The paper uses
thresholds of 1e-8 (start torus) and 1e-7 (scans).

## 12. Breakdown of tori and the limits of the families

**Smoothness and Fourier decay.** If K is analytic in a strip $|\mathrm{Im}\,\theta| < \delta$,
its Fourier coefficients decay like $e^{-2\pi\delta|k|}$. As a torus approaches breakdown, δ → 0:
the torus develops ever finer folds, and its coefficients decay ever more slowly. The number of
modes needed for a fixed accuracy grows roughly like $1/\delta$, without bound. At breakdown the
torus is no longer smooth (it may survive as a Cantor set, a "cantorus") and the parameterization
method can no longer represent it.

**The two signatures in the paper are two views of this:**

* *Orbit side:* the frequency drift of orbits rises, and the orbits stop being quasi-periodic
  (paper §6.2).
* *Torus side:* the highest Fourier mode needed rises sharply and anisotropically (paper §6.3).

Both happen at the same distance from the NRHO, about 12–20 km. That is the evidence that the
families end because the tori are breaking down, and not just because the solver gives up. The
computational response (anisotropic refinement, a relaxed floor) gained less than 1 km, which is
consistent with a hard physical edge.

**The surface is a different edge.** The CR3BP treats Enceladus as a point mass, so tori whose
QPOs dip below 252.1 km are still invariant tori of the model, but not usable orbits. For
low-perilune members (Id 97, Id 101) the surface ends the useful family well inside the regular
region.

## 13. How the other methods in the literature work

All of them solve $F(K(\theta)) = K(\theta+\omega)$ (or $\partial_\omega K = X(K)$ for a flow).
They differ in how K is represented, how the equation is solved, and whether the symplectic
structure is used.

| method | representation | solve | structure used | cost per step |
|---|---|---|---|---|
| Normal forms (Jorba 1999) | power series around the equilibrium or orbit | Lie-series transformations order by order | Hamiltonian | one-off; local, asymptotic |
| Fourier/Newton (Gómez–Mondelo 2001, Jorba 2001) | Fourier coefficients of an invariant curve | full Newton, dense Jacobian | none | $O(N^3)$ |
| Collocation/shooting (Olikara; Baresi–Olikara–Scheeres; McCarthy–Howell) | N points on an invariant circle of the stroboscopic map, propagated over T; rotation by DFT | Newton with phase conditions; pseudo-arclength continuation | none | $O(N^3)$ |
| Invariance PDE (Schilder–Osinga–Vogt 2005) | Fourier–Galerkin or finite differences for $\omega\cdot\partial_\theta u = f(u)$ | Newton plus pseudo-arclength | none (works for dissipative tori too) | $O(N^3)$ |
| **Parameterization method** (GJLV 2005; Haro–de la Llave 2006; Haro et al. 2016; Figueras–Haro–Luque 2017) | grid values and Fourier coefficients of K | quasi-Newton in the symplectic frame, cohomological equations | symplectic, Lagrangian | $O(N\log N)$ + N map calls |

Notes on each:

* **Normal forms** give the whole neighbourhood at once, with approximate first integrals, but
  the series are asymptotic. They have no error control on a given torus, and accuracy drops with
  distance.
* **Gómez–Mondelo** used Fourier/Newton with continuation to map out the tori around the collinear
  points. **Jorba 2001** also computes the linear normal behaviour (reducibility) of the curve.
* **Collocation** is the standard in astrodynamics, and it is what the 2023 QPO work used. It is
  robust and general. Olikara combines Gauss–Legendre collocation in time with Fourier in angle.
  Baresi et al. extend it to multiple shooting and families. McCarthy–Howell apply it to
  cislunar design.
* **Extensions of the parameterization method:** to flows, via the time-T flow map, handling the
  extra pair of degenerate directions from autonomy (Haro–Mondelo 2021); to quasi-periodically
  forced Hamiltonians (Fernández-Mora–Haro–Mondelo 2024) and whiskered tori with their
  convergence theory (2026); to partially integrable systems (Haro–Luque 2019); and to whiskered
  invariant circles in periodically perturbed planar CR3BPs (Kumar–Anderson–de la Llave 2022). For
  whiskered tori the normal bundle splits into centre, stable and unstable parts. The centre part
  is handled as above, and the hyperbolic parts are solved by contraction, with no small divisors.

## 14. Glossary

* **Action–angle variables.** Coordinates (I, φ) in which an integrable map is
  $(I,\varphi)\mapsto(I,\varphi+\omega(I))$.
* **Cantor family.** A set with no intervals but positive measure. The KAM tori are indexed by
  one.
* **Cohomological equation.** $u(\theta) - u(\theta+\omega) = v(\theta)$, diagonal in Fourier
  space with divisors $1 - e^{2\pi ik\cdot\omega}$.
* **Diophantine.** $|k\cdot\omega - m| \ge \gamma|k|^{-\tau}$: badly approximable by rationals.
* **Elliptic (fixed point or pair).** Eigenvalues on the unit circle, so the linearization is a
  rotation.
* **Exact symplectic.** Symplectic with $F^*\alpha - \alpha$ exact, where $\Omega = d\alpha$.
  This rules out tori that "drift".
* **Frequency drift.** The difference between rotation numbers from two halves of an orbit, used
  as a regularity indicator.
* **Invariance error.** $E = F\circ K - K(\cdot+\omega)$, reported as a sup-norm over the grid
  in LU.
* **Lagrangian torus.** $DK^\top\Omega DK = 0$: a half-dimensional torus with zero symplectic
  area on every tangent plane.
* **Linear rotation numbers ρ_lin.** The rotation angles over 2π of the elliptic eigenvalues of
  DF at the NRHO.
* **Low-pass filter.** Zeroing the untrustworthy high-|k| half of the spectrum after each step.
* **Mode mix (c₁, c₂).** The relative amplitudes of the two elliptic modes in the start torus.
  It selects the path through the two-parameter family.
* **NAFF.** Numerical Analysis of Fundamental Frequencies (Laskar): windowed Fourier refinement
  of the frequencies of a quasi-periodic signal.
* **Noise floor.** The accuracy of F itself, below which Newton cannot reduce the error.
* **Rotation vector ω.** The translation in the angles per map iterate on the torus, in cycles.
* **Symplectic map.** $DF^\top\Omega DF = \Omega$: preserves the symplectic 2-form.
* **Torsion T(θ).** The upper-right block of the linearized map in the (L, N) frame. Its average
  ⟨T⟩ must be invertible (the twist condition).
* **Twist condition.** $\det\langle T\rangle\neq 0$, or, in the integrable picture,
  $\det(\partial\omega/\partial I)\neq 0$.
