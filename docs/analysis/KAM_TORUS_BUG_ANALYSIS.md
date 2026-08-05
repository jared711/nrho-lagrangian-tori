# KAM Torus Bug: Comprehensive Line-by-Line Analysis

## Executive Summary

**Bug Location**: param.cc lines 264-266 (duplicate kam_torus calls)  
**Root Causes Identified**: 
1. Duplicate kam_torus() calls compound Fourier truncation errors
2. Missing angle offset in map evaluation (commented out in commit ca0df7c, line 613-614)
3. In-place modification of paramR/paramF without state preservation

**Impact**: Second kam_torus() call INCREASES error instead of decreasing it, defeating the KAM correction algorithm.

---

## Background: The Parameterization Method

Your code implements **Haro & Luque's parameterization method** (from "The parameterization method for invariant manifolds," 2016):

- **Goal**: Compute an invariant torus K(θ) in phase space where θ = (θ₁, θ₂) are angle variables
- **Invariance Equation**: K(θ + ω) = Map(K(θ))
  - If K satisfies this, following the map forward is equivalent to rotating angles by ω
- **Method**: Represent K as Fourier series, use quasi-Newton to minimize invariance error
- **Correction**: K_new = K + L·ξ_L + N·ξ_N
  - L = tangent frame (derivatives ∂K/∂θ)
  - N = normal frame (symplectic complement)
  - ξ_L, ξ_N = corrections solving cohomological equations

---

## The Calling Code (Lines 264-266)

```cpp
263 /* If the initial torus is not invariant, we should uncomment the following line: */
264 conv = kam_torus(paramR,paramF,omega,error,nn,nelem,tail0,tails,2,
                     map_CR3BP,sform_CR3BP,gform_CR3BP,normal0_CR3BP);
265
266 conv = kam_torus(paramR,paramF,omega,error,nn,nelem,tail0,tails,2,
                     map_CR3BP,sform_CR3BP,gform_CR3BP,normal0_CR3BP);
```

**Problem**: Two identical calls with no state check between them.

**Arguments**:
- `paramR`, `paramF`: Torus parameterization (real-space and Fourier), **passed by reference** (modified in-place)
- `omega`: Rotation frequencies
- `error`: L1 norm of invariance error (output)
- `nn`: Grid dimensions
- `nelem`: Total grid points
- `tail0`, `tails`: Fourier tail convergence flags
- `Case=2`: Use metric-based normal frame (lines 709-715)
- Function pointers: CR3BP map, symplectic form, metric, normal0

---

## kam_torus() Line-by-Line Analysis

### Function Signature (Lines 492-497)

```cpp
int kam_torus(matrix &paramR, matrix &paramF, myreal *omega, myreal &error, int *nn, int nelem,
              int &tail0, int *tails, int Case,
              void (*map)(complex *, complex *, complex **, complex *),
              void (*sform)(complex *, complex **),
              void (*gform)(complex *, complex **),
              void (*normal0)(matrix &, int *, int))
```

**Key observation**: `paramR` and `paramF` are **references** - modifications persist after return.

**Return values**:
- `0`: Not converged, needs more iterations
- `1`: Converged (error < tolinva), line 685

---

### STEP 0: Tail Evaluation (Lines 575-603)

```cpp
575-578  tail(paramF, tg); // Compute Fourier tail size (eq 4.104 in Haro et al. 2016)
580      tail0 = 0;
588-597  for (int i = 0; i < DTOR; i++) {
             if (tg[i] > toltail) {  // If tail too large, stop
                 tails[i] = 1;
                 tail0 = 1;
             }
         }
599-600  if (tail0 == 1)
             return 0;  // Grid resolution insufficient, need more Fourier modes
601      clean(paramF);  // ⚠️ CRITICAL: Truncates high-frequency modes
602      paramR = fft_B(paramF);  // ⚠️ OVERWRITES input with cleaned version
```

**⚠️ BUG MECHANISM #1: Compounding Truncation**

Line 601-602 are executed **every time** kam_torus() is called:
1. **First call (line 264)**: 
   - `clean(paramF)` removes small Fourier coefficients
   - `paramR = fft_B(paramF)` overwrites the input
   - Correction is computed and applied (line 823-824)
   - Returns with **truncated** paramR/paramF

2. **Second call (line 266)**:
   - Receives **already-truncated** paramR/paramF from first call
   - `clean(paramF)` removes MORE small coefficients
   - Correction is computed on **degraded data**
   - **Accumulated truncation error** dominates, correction goes wrong direction

**Why this matters**: Small Fourier coefficients encode fine structure needed to compute the correct tangent/normal frames (L, N). Truncating twice loses information critical for the Newton step.

---

### STEP 1: Invariance Error (Lines 604-651)

```cpp
607-615  for (int l = 0; l < nelem; l++) {  // Loop over all grid points
             indices(l, nn, index, DTOR);  // Convert linear index l to multi-index (k₁, k₂)
             for (int i = 0; i < DMAP; i++) {
                 z[i] = paramR.coef[i][0].elem[l];  // Extract point K(θ_k)
                 // ⚠️ BUG #2: These lines were commented out in commit ca0df7c:
                 // if (i < DTOR)
                 //     z[i] = z[i] + ((double)index[i]) / ((double)nn[i]);
             }
616          (*map)(z, fz, Dfz, depfz);  // Compute Map(K(θ_k)) and its derivative
617-619      (*sform)(z, Omegaz);  // Symplectic form
             if (Case != 1)
                 (*gform)(z, Metricz);  // Riemannian metric
620-630      // Store Map(K), DMap, Ω, g in grid format
         }
```

**⚠️ BUG MECHANISM #2: Missing Angle Offset**

The commented-out lines (613-614) added `index[i]/nn[i]` to the position coordinates before map evaluation.

**What this does**:
- `index[i]` = grid index (0, 1, ..., nn[i]-1)
- `index[i]/nn[i]` = normalized angle coordinate (0, 1/n, 2/n, ..., (n-1)/n)
- For Fourier discretization: θ_k = 2π(k₁/n₁, k₂/n₂)

**Why it's needed**:
- **Hypothesis**: The parameterization K(θ) might be stored as *relative* coordinates
- The map evaluation needs *absolute* position = K(θ) + θ
- Without this offset, the map is evaluated at the **wrong points**
- This causes the invariance error E = Map(K) - K∘shift to be computed incorrectly
- The correction then points in the **wrong direction**, increasing error instead of decreasing it

**Evidence from git history**:
- Commit ca0df7c message: "Jared wants to comment out these lines on 6/27/24"
- Same commit message: "torus correction seems to be wrong"
- Correlation suggests this change introduced the bug

---

```cpp
633-643  KshiftF = shift(paramF, omega);  // Shift Fourier modes: K̂(k) → K̂(k)·e^{2πi k·ω}
         KshiftR = fft_B(KshiftF);       // Convert back to real space
         
         // ⚠️ Note: Angle offset IS applied here (line 641):
         for (int l = 0; l < nelem; l++) {
             indices(l, nn, index, DTOR);
             for (int i = 0; i < DTOR; i++) {
                 KshiftR.coef[i][0].elem[l] = KshiftR.coef[i][0].elem[l] 
                     + ((double)index[i]) / ((double)nn[i]) + omega[i];
             }
         }

645      ErrorR = FparamR - KshiftR;  // Invariance error in real space
647-648  ErrorF = fft_F(ErrorR);      // Convert to Fourier space
         error = norm(ErrorF);         // L1 norm: sum |ErrorF_k|
```

**Key asymmetry**: 
- Angle offset IS added to KshiftR (line 641)
- Angle offset NOT added to z before map evaluation (line 613-614, commented out)
- This asymmetry causes E = Map(K(θ)) - K(θ+ω) to have a systematic bias

---

```cpp
652-686  if (error < tolinva) {  // If already converged...
             cout << "#     - No correction is needed!" << endl;
             // Clean up memory
             return 1;  // ⚠️ Only place where return 1 (converged)
         }
```

**Observation**: No convergence check between the duplicate calls on lines 264-266.
- If first call reduces error below `tolinva`, it returns 1
- But the calling code doesn't check `conv` before line 266
- Second call will still run even if first converged

---

### STEP 2: Symplectic Frame Construction (Lines 688-724)

```cpp
692-693  DparamF = diff(paramF);      // Differentiate: ∂K/∂θ in Fourier space
         DparamR = fft_B(DparamF);   // Convert to real space

695      LR = DparamR;  // Tangent frame L = DK (derivatives along torus)

709-715  else if (Case == 2) {  // Metric-based normal (used in your calls)
             GR = trans(LR) * MetricKR * LR;  // Gram matrix G = L^T g L
             BR = inv(GR);                     // B = G^{-1}
             AR = trans(BR) * trans(LR) * MetricKR * inv(OmegaKR) * MetricKR * LR * BR * val05;
             NR = LR * AR - inv(OmegaKR) * MetricKR * LR * BR;  // Normal frame
         }

722-723  LF = fft_F(LR);  // Convert L to Fourier
         NF = fft_F(NR);  // Convert N to Fourier
```

**Purpose**: Build tangent-normal frame (L, N) satisfying:
- L^T Ω N = 0 (symplectic orthogonality)
- L, N span phase space at each point

**Numerical issues**:
- Line 712: `inv(GR)` - if G is near-singular, large errors amplify here
- Line 713: `inv(OmegaKR)` - symplectic form should be non-singular, but numerical errors can accumulate
- If the cleaned paramF lost important modes, L and N will be inaccurate

---

### STEP 3: Cohomological Equations (Lines 725-788)

```cpp
728-731  LshiftF = shift(LF, omega);  // Shift L and N
         LshiftR = fft_B(LshiftF);
         NshiftF = shift(NF, omega);
         NshiftR = fft_B(NshiftF);

732-734  OmegaKF = fft_F(OmegaKR);
         OmegaKshiftF = shift(OmegaKF, omega);
         OmegaKshiftR = fft_B(OmegaKshiftF);

735-736  etaLR = -trans(NshiftR) * OmegaKshiftR * ErrorR;  // Project error onto L
         etaNR = trans(LshiftR) * OmegaKshiftR * ErrorR;   // Project error onto N

737      twistR = trans(NshiftR) * OmegaKshiftR * DFKR * NR;  // "Twist" matrix
```

**Cohomological equation for normal component**:
```cpp
738-740  etaNF = fft_F(etaNR);
         RetaNF = cohomological(etaNF, omega);  // Solve: ω·∇ξ_N = η_N (in Fourier)
         RetaNR = fft_B(RetaNF);
```

The `cohomological()` function inverts the operator ω·∇ in Fourier space:
- ξ_N(k) = η_N(k) / (2πi k·ω)  (small divisor problem!)
- When k·ω ≈ 0 (resonance), division amplifies errors
- For typical NRHOs, ω is irrational → no exact resonances, but near-resonances exist

**Twist matrix inversion**:
```cpp
741      newetaR = etaLR - twistR * RetaNR;  // Modified equation for tangent

743-758  aver(twistR, twist0);  // Average twist over grid → DTOR×DTOR matrix
         tolqr = tolinte;
         qrdcmp(twist0, DTOR, DTOR, tolqr);  // QR decomposition for inversion
         // ... solve twist0 * xiN0 = neweta0 ...
         for (int i = 0; i < DTOR; i++) {
             invT[i][j] = solc[i];  // Store inverse twist
             global_twist = global_twist + invT[i][j];  // Condition indicator
         }
```

**⚠️ POTENTIAL INSTABILITY**:
- Line 745: `qrdcmp()` - if twist matrix is near-singular, QR decomposition is unstable
- Line 756: `global_twist` sums all entries of inv(twist) - large value = ill-conditioned
- If truncation from `clean()` degraded the data, twist matrix can become singular

```cpp
774-782  xiNR = RetaNR;
         for (int l = 0; l < nelem; l++) {
             for (int i = 0; i < DTOR; i++) {
                 xiNR.coef[i][0].elem[l] = xiNR.coef[i][0].elem[l] + xiN0[i][0];
             }
         }  // Add constant part to periodic part

784-787  newetaR = etaLR - twistR * xiNR;  // Update residual
         newetaF = fft_F(newetaR);
         xiLF = cohomological(newetaF, omega);  // Solve for tangent correction
         xiLR = fft_B(xiLF);
```

---

### STEP 4: Apply Correction (Lines 789-834)

```cpp
792-793  newparamR = paramR + LR * xiLR + NR * xiNR;  // K_new = K + L·ξ_L + N·ξ_N
         newparamF = fft_F(newparamR);

823-824  paramR = newparamR;  // ⚠️ OVERWRITE input references
         paramF = newparamF;

826-831  /* Just to show the size of the correction */
         newparamR = LR * xiLR + NR * xiNR;
         newparamF = fft_F(newparamR);
         aux = norm(newparamF);
         cout << "#     - Norm of the correction: " << aux << endl;

833      return 0;  // Not converged, need more iterations
```

**⚠️ BUG MECHANISM #3: In-Place Modification**

Lines 823-824 overwrite the reference parameters. When line 266 calls kam_torus() again:
- Input is the **already-corrected** torus from line 264
- Line 601-602 **truncate it again**
- Corrections compound, but so do truncation errors
- Net effect: error increases instead of decreases

---

## Root Cause Summary

### Primary Bug: Duplicate Call with Truncation
Lines 264-266 call kam_torus() twice on the same data:
1. First call cleans/truncates Fourier modes (line 601-602)
2. Corrections are applied (line 823-824)
3. Second call cleans/truncates **again** on already-degraded data
4. Accumulated truncation error dominates → divergence

**Fix**: Remove line 266 (duplicate call).

---

### Secondary Bug: Missing Angle Offset
Lines 613-614 (commented out) should add `index[i]/nn[i]` to coordinates:
- Creates asymmetry in error computation
- Map evaluated at wrong points
- Correction points in wrong direction

**Fix Options**:
1. **Restore the angle offset** (uncomment lines 613-614) - but this contradicts the comment "Jared wants to comment out"
2. **Investigate why it was commented out** - was there a reason?
3. **Test both versions** - run with and without offset, compare convergence

---

### Tertiary Issue: No Convergence Check
Line 266 calls kam_torus() unconditionally:
- Doesn't check if line 264 converged (`conv == 1`)
- Doesn't check if error increased
- No adaptive step size or damping

**Fix**: Add convergence logic:
```cpp
conv = kam_torus(paramR,paramF,omega,error,nn,nelem,tail0,tails,2,
                 map_CR3BP,sform_CR3BP,gform_CR3BP,normal0_CR3BP);
if (conv == 0 && error > some_threshold) {  // If not converged and error reasonable...
    conv = kam_torus(paramR,paramF,omega,error,nn,nelem,tail0,tails,2,
                     map_CR3BP,sform_CR3BP,gform_CR3BP,normal0_CR3BP);
}
```

---

## Recommended Debugging Steps

### 1. Immediate Fix (High Confidence)
**Remove line 266** - eliminate the duplicate call:
```cpp
// Before:
264 conv = kam_torus(...);
265
266 conv = kam_torus(...);  // ← DELETE THIS LINE

// After:
264 conv = kam_torus(...);
```

**Rationale**: Compounding truncation is the most obvious bug. Single call should work if initial guess is good.

---

### 2. Test Angle Offset (Medium Confidence)
**Restore lines 613-614** and test:
```cpp
for (int i = 0; i < DMAP; i++) {
    z[i] = paramR.coef[i][0].elem[l];
    if (i < DTOR)  // ← UNCOMMENT
        z[i] = z[i] + ((double)index[i]) / ((double)nn[i]);  // ← UNCOMMENT
}
```

**Test cases**:
- Run with offset: `./param approxQPO.csv` (original version before ca0df7c)
- Run without offset: current version
- Compare: error trajectory, convergence rate, final accuracy

**Look for**:
- Does error decrease monotonically with offset?
- Does final torus have better invariance error?

---

### 3. Instrumentation (Understand Failure Mode)
Add diagnostic output to kam_torus() to track:
```cpp
// After line 648:
cout << "#     - Error of invariance: " << error << endl;
cout << "#       DEBUG: Norm of paramF before clean: " << norm(paramF) << endl;

// After line 601:
clean(paramF);
cout << "#       DEBUG: Norm of paramF after clean: " << norm(paramF) << endl;

// After line 761:
cout << "#     - Norm inverse twist: " << global_twist << endl;
cout << "#       DEBUG: Is twist ill-conditioned? " << (global_twist > 1e6 ? "YES" : "no") << endl;

// After line 831:
cout << "#     - Norm of the correction: " << aux << endl;
cout << "#       DEBUG: Correction size relative to torus: " << (aux / norm(paramF)) << endl;
```

**Look for**:
- Is `clean()` removing a significant fraction of the norm?
- Is `global_twist` growing large (> 10^6) indicating near-singularity?
- Is correction size comparable to torus size (bad) or small (good)?

---

### 4. Continuation Loop Investigation
Lines 271-273 save the corrected torus:
```cpp
271 paramR0 = matrix(paramR);
272 paramF0 = matrix(paramF);
273 epsilon0 = epsilon;
```

**Question**: Where is the continuation loop that actually uses kam_torus() iteratively?
- The duplicate call (lines 264-266) suggests someone tried to "iterate twice"
- But proper continuation should be a loop with convergence checks

**Search for**: A loop structure like:
```cpp
while (error > tolerance && iter < maxiter) {
    conv = kam_torus(...);
    if (conv == 1) break;  // Converged
    iter++;
}
```

If this doesn't exist, it needs to be implemented.

---

## References

### Papers (from Alex Haro literature search)
1. **Haro & de la Llave (2006)**: "A parameterization method for the computation of invariant tori and their whiskers in quasi-periodic maps: numerical algorithms"
   - Foundational paper for this method

2. **Haro (2016)**: "The Parameterization Method for Quasi-Periodic Systems: From Rigorous Results to Validated Numerics"
   - Section 4.104: Tail definition
   - Algorithm 4.32: KAM correction procedure

3. **Haro, Fernández-Mora, Mondelo (2024)**: "On the convergence of flow map parameterization methods in Hamiltonian systems"
   - When/why Newton schemes fail - relevant for debugging divergence

### Code Comments
- Line 578: References eq 4.104 in Haro et al. (2016)
- Line 263: "If the initial torus is not invariant, we should uncomment the following line"
  - Suggests lines 264-266 were meant to be optional, not both active

### Git History
- Commit ca0df7c: "Poincaré map works, but torus correction is wrong"
  - Commented out angle offset (lines 613-614)
  - Bug correlates with this change

---

## Questions for User

1. **Why were lines 613-614 commented out?** (commit ca0df7c, "Jared wants to comment out these lines")
   - Was there a specific problem it caused?
   - Did it fix something else but break torus correction?

2. **What is the expected convergence rate?**
   - How many kam_torus() calls should it take to converge?
   - What's a typical initial error vs final error?

3. **Is there supposed to be a continuation loop?**
   - Lines 264-266 suggest manual "two iterations"
   - Should there be a while loop instead?

4. **What are the tolerance values?**
   - `toltail`: Fourier tail threshold
   - `tolinva`: Invariance error threshold for convergence
   - `tolinte`: QR decomposition tolerance
   - Are these set appropriately for your problem?

---

## Next Steps

**Immediate**:
1. Remove line 266 (duplicate call)
2. Run `make clean && make`
3. Test: `./param approxQPO.csv`
4. Check if single kam_torus() call reduces error

**Short-term**:
1. Decide on angle offset (lines 613-614) - test with/without
2. Add diagnostic output to track error, truncation, conditioning
3. Verify convergence criteria and tolerance values

**Medium-term**:
1. Implement proper continuation loop with adaptive iterations
2. Add error handling for ill-conditioned twist matrix
3. Consider adaptive Fourier mode selection (auto-refine grid if tail too large)

**Long-term** (for publication):
1. A posteriori validation (Haro 2016 method)
2. Compare results to reference implementations (if any found)
3. Document algorithm choices and numerical parameters
