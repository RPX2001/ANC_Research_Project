# ANC Simulation Error Summary and Fixes

## Errors Found and Fixed

### 1. **Critical: Variable Definition Order in anc_system.m**
**Problem:** The script tried to use variables before they were defined.
- Lines 13-30 attempted to load audio files using `Ns_target` and `fs` which were not yet defined
- Lines 56+ tried to use `h_pm`, `h_sm`, `h_pee`, `h_see`, `retm_est`, `nfft`, `win`, `wlen`, `hop` which were never defined in the file

**Fix:** Restructured `anc_system.m` to:
- Run AFTER `new_anc_simulation.m` (which provides all required variables)
- Added validation to check if required variables exist in workspace
- Removed redundant source signal loading (already done in `new_anc_simulation.m`)
- Proper parameters defined before use

### 2. **Critical: Wrong Loop Variable Indexing in new_anc_simulation.m**  
**Problem:** Three FFT computation loops used the wrong index variable:

**Line ~224:** 
```matlab
for k = 1:num_pri_src
    Ht = squeeze(h_pe(:,:,l)).';    % BUG: 'l' doesn't exist, should be 'k'
```

**Line ~234:**
```matlab
for l = 1:num_sec_src
    Ht = squeeze(h_sm(:,:,k)).';    % BUG: should be 'l' not 'k'
```

**Line ~244:**
```matlab
for l = 1:num_sec_src
    Ht = squeeze(h_se(:,:,k)).';    % BUG: should be 'l' not 'k'
```

**Impact:** These bugs would cause:
- Incorrect frequency responses for paths
- Wrong ReTM estimation  
- Failed ANC performance
- Possibly no PSD plots or incorrect results

**Fix:** Changed to use correct loop variables (`k` for primary sources, `l` for secondary sources)

### 3. **Missing Integration Between Scripts**
**Problem:** `anc_system.m` was written as standalone but required outputs from `new_anc_simulation.m`

**Fix:** Added `run('anc_system.m')` at the end of `new_anc_simulation.m` to create complete workflow

### 4. **Code Organization Issues**
**Problem:** 
- Poor documentation
- No progress reporting
- Inconsistent variable naming
- No validation of intermediate results

**Fix:** Added:
- Clear section headers
- Progress reporting during ANC loop  
- Variable existence validation
- Performance summary metrics
- Better plot labels and titles

---

## How to Run the Corrected Code

### Step 1: Navigate to the simulation directory
```matlab
cd('d:\Work\ANC_Research_Project\Simulation_ANC')
```

### Step 2: Run the complete simulation
```matlab
clear all; close all;
run('new_anc_simulation.m');
```

This will:
1. Generate room impulse responses (RIRs)
2. Load and prepare source signals
3. Generate noise signals  
4. Compute monitoring and error mic signals
5. Estimate the Relative Transfer Matrix (ReTM)
6. Verify ReTM estimation accuracy
7. **Automatically run ANC control stage**
8. Generate three plots:
   - Time-domain comparison (ANC OFF vs ON)
   - Power Spectral Density comparison
   - Noise Reduction vs Frequency
9. Display performance metrics

### Step 3: Check the results
You should now see:
- **Figure 1-2:** ReTM verification plots (from tuning stage)
- **Figure 3:** Time-domain comparison at evaluation mic
- **Figure 4:** PSD comparison (dB/Hz) - **THIS WAS MISSING BEFORE**
- **Figure 5:** Noise reduction vs frequency (dB) - **THIS MAY HAVE HAD ERRORS BEFORE**

### Expected Console Output:
```
======================================== 
   Tuning Stage Complete
   Starting ANC Control Stage
========================================

=== Starting ANC Control Stage ===
System dimensions: K=2 sources, M=8 mon mics, L=2 speakers, Neval=2 eval mics
ReTM extracted: 4 virtual mics mapped from 8 monitoring mics

=== Initializing ANC Controller ===
  Step size mu = 3.0e-04
  Adapting 1025 frequency bins (0.0 - 4000.0 Hz)

=== Running ANC Control Loop ===
  Processing frames 51 to ...
  Progress: 0% (frame 51/...)
  ...
  Progress: 100% (frame .../...)
  ANC loop completed in XX.X seconds

=== Performance Summary ===
  Average noise reduction (50-600 Hz): XX.X dB
  Peak noise reduction: XX.X dB at XXX.X Hz

=== ANC Analysis Complete ===
```

---

## Understanding the Results

### Power Spectral Density (PSD) Plot
- **Blue line:** Noise spectrum BEFORE ANC (baseline)
- **Red line:** Noise spectrum AFTER ANC  
- **Expected:** Red line should be below blue line in the control frequency band (50-600 Hz)

### Noise Reduction Plot  
- Shows frequency-dependent performance
- **Positive values:** ANC is reducing noise at that frequency
- **Negative values:** ANC is amplifying noise (should be minimal)
- **NaN gaps:** Frequencies with very low input energy (ignored)

---

## Troubleshooting

### If PSD plot still doesn't appear:
1. Check that `eval_off` and `eval_on` contain non-zero data
2. Verify `fs = 8000` is correct
3. Check MATLAB's figure windows aren't hidden

### If noise reduction is poor:
1. Try adjusting `mu` (learning rate) in anc_system.m:
   - Increase for faster convergence: `mu = 5e-4`
   - Decrease for stability: `mu = 1e-4`
2. Check the frequency adaptation range
3. Verify ReTM estimation is accurate (check earlier plots)

### If errors occur about missing variables:
- Make sure to run `new_anc_simulation.m` completely, not `anc_system.m` directly
- Don't clear workspace between scripts

---

## Files Modified
1. **anc_system.m** - Complete restructure
2. **new_anc_simulation.m** - Fixed indexing bugs, added integration

## Files Created
1. **ERROR_SUMMARY.md** - This file

---

## Next Steps for Further Improvement
1. Tune `mu` (step size) for optimal convergence speed vs stability
2. Experiment with different adaptation frequency ranges
3. Try different control time windows (`t_start`, `t_end`)
4. Compare results with paper figures  
5. Test with different source signals
6. Analyze convergence behavior over time

---

**Date:** February 12, 2026  
**Status:** ✅ All critical errors fixed, code ready to run
