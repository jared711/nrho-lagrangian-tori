# Continuation Loop Fix

**Date**: 2026-08-05  
**Issue**: Continuation loop runs indefinitely without termination criteria  
**Fix**: Added maximum epsilon bound and maximum step counter  

---

## Problem

The continuation loop in `param.cc` (lines 330-432) had no upper bounds:

```cpp
do {
    epsilon = epsilon0 + deps;  // Increment epsilon
    
    // Newton iterations to compute torus...
    
    if (conv == 1) {
        epsilon0 = epsilon;  // Success, save and continue
    } else {
        deps = deps / 10;    // Failure, reduce step size
        if (deps < 1e-5)
            fail = 1;        // Only exit condition!
    }
} while (fail == 0);
```

**Issues**:
1. Loop only exits if Newton method fails AND deps < 1e-5
2. No upper bound on epsilon
3. No limit on number of continuation steps
4. For successful continuation, would run forever
5. Each continuation step takes ~2 hours for 128×128 grid

**Observed Behavior**:
- Started at epsilon = 0
- First kam_torus() call (line 264, before continuation) completed successfully
- Entered continuation loop, started computing epsilon = 5e-5
- After 2 hours, still in first continuation step computing Poincaré maps
- Manually stopped after 118 minutes

---

## Solution

Added two termination criteria:

### 1. Maximum Epsilon Bound

```cpp
const double MAX_EPSILON = 0.01;  // Stop when epsilon exceeds this

if (epsilon > MAX_EPSILON) {
    cout << "# Reached maximum epsilon = " << epsilon << " (limit: " << MAX_EPSILON << ")" << endl;
    cout << "# Terminating continuation." << endl;
    fail = 1;
    break;
}
```

**Rationale**: 
- Epsilon represents perturbation parameter in continuation
- Starting at epsilon=0 (nominal NRHO), continuing to epsilon=0.01 provides reasonable coverage
- Value is adjustable based on application needs

### 2. Maximum Continuation Steps

```cpp
const int MAX_CONT_STEPS = 10;  // Maximum number of continuation steps

if (cont_step >= MAX_CONT_STEPS) {
    cout << "# Reached maximum continuation steps = " << cont_step << " (limit: " << MAX_CONT_STEPS << ")" << endl;
    cout << "# Terminating continuation." << endl;
    fail = 1;
    break;
}
```

**Rationale**:
- Each step takes hours for high-resolution grids
- 10 steps provides reasonable parameter coverage
- Prevents runaway continuation even if epsilon bound not reached

### 3. Improved Progress Reporting

```cpp
cont_step++;
cout << "# Continuation step " << cont_step << " / " << MAX_CONT_STEPS << endl;
cout << "# We try to compute the torus for epsilon=" << epsilon << " (limit: " << MAX_EPSILON << ")" << endl;
```

**Benefits**:
- User can see progress toward termination
- Clear indication of how many steps remain
- Shows current epsilon value relative to limit

---

## Modified Code

**Location**: `param.cc` lines 329-360

### Before:
```cpp
/**** START Continuation with respect to epsilon ****/
int fail = 0;
do
{
    paramR = paramR0;
    paramF = paramF0;
    epsilon = epsilon0 + deps;

    /**** START Newton method to correct the invariant torus ****/
    cout << "# We try to compute the torus for epsilon=" << epsilon << endl;
    iter = 0;
    do
    {
        ...
    } while (conv == 0 && iter < MNEW && tail0 == 0);
    ...
} while (fail == 0);
```

### After:
```cpp
/**** START Continuation with respect to epsilon ****/
int fail = 0;
int cont_step = 0;
const int MAX_CONT_STEPS = 10;  // Maximum number of continuation steps
const double MAX_EPSILON = 0.01; // Maximum value of epsilon to continue
do
{
    paramR = paramR0;
    paramF = paramF0;
    epsilon = epsilon0 + deps;

    // Check continuation termination criteria
    if (epsilon > MAX_EPSILON) {
        cout << "# Reached maximum epsilon = " << epsilon << " (limit: " << MAX_EPSILON << ")" << endl;
        cout << "# Terminating continuation." << endl;
        fail = 1;
        break;
    }
    if (cont_step >= MAX_CONT_STEPS) {
        cout << "# Reached maximum continuation steps = " << cont_step << " (limit: " << MAX_CONT_STEPS << ")" << endl;
        cout << "# Terminating continuation." << endl;
        fail = 1;
        break;
    }
    cont_step++;

    /**** START Newton method to correct the invariant torus ****/
    cout << "# Continuation step " << cont_step << " / " << MAX_CONT_STEPS << endl;
    cout << "# We try to compute the torus for epsilon=" << epsilon << " (limit: " << MAX_EPSILON << ")" << endl;
    iter = 0;
    do
    {
        ...
    } while (conv == 0 && iter < MNEW && tail0 == 0);
    ...
} while (fail == 0);
```

---

## Customization

Users can adjust these parameters based on their needs:

### For Quick Testing:
```cpp
const int MAX_CONT_STEPS = 1;      // Just one continuation step
const double MAX_EPSILON = 0.0001; // Very small perturbation
```

### For Thorough Exploration:
```cpp
const int MAX_CONT_STEPS = 50;     // Many continuation steps
const double MAX_EPSILON = 0.1;    // Larger perturbation range
```

### For Single Torus Computation (No Continuation):
```cpp
const int MAX_CONT_STEPS = 0;      // Skip continuation entirely
```

Or comment out the entire continuation loop and use only the initial kam_torus() call on line 264.

---

## Expected Runtime

With these limits:

**Scenario 1**: All steps converge successfully
- Runtime: ~MAX_CONT_STEPS × (time per step)
- For 128×128 grid: 10 steps × 2 hours = ~20 hours

**Scenario 2**: Hit epsilon limit early
- Runtime: (epsilon_limit / deps) × (time per step)
- For MAX_EPSILON=0.01, deps=5e-5: 0.01/5e-5 = 200 steps
- But MAX_CONT_STEPS=10 will terminate first → ~20 hours

**Scenario 3**: Newton method fails
- Runtime: Until failure, then deps reduced and retry
- Will eventually hit deps < 1e-5 or step limit

**Recommendation**: For initial testing, use MAX_CONT_STEPS=1 to verify continuation works without waiting hours.

---

## Testing

### Quick Test (Recommended First):
```cpp
const int MAX_CONT_STEPS = 1;
const double MAX_EPSILON = 0.01;
```

Run:
```bash
make
./param approxQPO.csv 2>&1 | tee param_test_cont1.txt
```

Expected runtime: ~2-4 hours (one continuation step after initial torus)

### Full Test:
```cpp
const int MAX_CONT_STEPS = 10;
const double MAX_EPSILON = 0.01;
```

Expected runtime: ~20 hours or until epsilon limit reached

---

## Status

✅ **Fix implemented**: Lines 329-360 in param.cc  
✅ **Code compiled**: param executable rebuilt  
⏳ **Testing**: Not yet tested with new limits  
📝 **Documentation**: This file + TEST_RESULTS.md updated  

---

## Next Steps

1. **Test with MAX_CONT_STEPS=1** to verify continuation works
2. **Monitor output** for proper termination messages
3. **Adjust limits** based on results and needs
4. **Commit changes** once verified working

---

## Related Files

- `param.cc` - Main implementation (continuation loop lines 329-432)
- `TEST_RESULTS.md` - Test results and bug verification
- `KAM_TORUS_BUG_ANALYSIS.md` - Original bug analysis
- `PAPER_SUMMARY.md` - Paper updated with results
