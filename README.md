# ANC Research Project

Active Noise Cancellation (ANC) simulation project focused on virtual sensing in enclosed spaces, with a comparison between RM and ReTM estimation methods.

## Project Overview

This repository contains:

- MATLAB simulation code for static-noise ANC scenarios.
- Room impulse response (RIR) generation utilities.
- ReTM/RM estimation functions.
- Scripts to evaluate ANC performance using PSD and noise reduction plots.

## Code Structure

### Top-Level

- `README.md`: this file.
- `documents/`: research notes, method documentation, and result artifacts.
- `Simulation_ANC/`: main simulation and algorithm implementation.

### `Simulation_ANC/`

- `new_anc_simulation.m`
	- Main end-to-end script for RM-based virtual sensing ANC.
	- Performs room/acoustic setup, tuning-stage signal generation, RM estimation, verification plots, and then runs control stage via `anc_system.m`.

- `anc_system.m`
	- ANC control stage (frequency-domain FxLMS with virtual sensing map from `rm_est`).
	- Produces time-domain, PSD, and noise-reduction outputs at evaluation microphones.

- `compare_rm_retm_psd.m`
	- Runs the same scenario twice (RM then ReTM), reusing the same ANC pipeline.
	- Plots combined RM vs ReTM PSD and noise-reduction curves.
	- Saves results to `rm_retm_comparison_results.mat`.

- `ERROR_SUMMARY.md`
	- Notes on previously identified issues and fixes.

- `ReTM/`
	- Core estimators and examples:
	- `rm_estimate.m`, `retm_estimate.m`, `retf_estimate.m`
	- `example_retm_denoise.m`, `example_retf_denoise.m`
	- `+shaasp/`: signal processing and acoustic helper functions.

- `RIR_Generator/`
	- `rir_generator.mexw64` plus examples and reference docs.
	- Used by simulation scripts to generate acoustic paths.

- `Sound_sources/`
	- Input audio used by the simulation (for example `factory1.wav`, `factory2.wav`).

## Requirements

- MATLAB (recommended recent version).
- Signal Processing Toolbox (used by `stft`, `istft`, `pwelch`, `bandpass`, `hann`).
- Windows-compatible prebuilt RIR MEX binary exists at:
	- `Simulation_ANC/RIR_Generator/rir_generator.mexw64`

If you are on another OS, you may need to compile `rir_generator.cpp` for your platform.

## How To Run The Algorithm

Open MATLAB, then run from the `Simulation_ANC` folder so relative paths resolve correctly.

```matlab
cd('d:/Work/ANC_Research_Project/Simulation_ANC')
clearvars; close all; clc;
```

### Option 1: Run Main ANC Pipeline (RM-based)

```matlab
run('new_anc_simulation.m')
```

What this does:

1. Generates room impulse responses.
2. Loads source audio and creates tuning/control signals.
3. Builds monitoring and error microphone signals.
4. Estimates RM (`rm_est`) from tuning data.
5. Runs ANC control stage through `anc_system.m`.
6. Produces plots and summary metrics:
	 - Evaluation mic time-domain comparison.
	 - PSD (ANC OFF vs ANC ON).
	 - Noise reduction vs frequency.

### Option 2: Compare RM vs ReTM Using Same Scenario

```matlab
run('compare_rm_retm_psd.m')
```

What this does:

1. Executes RM pipeline via `new_anc_simulation.m`.
2. Re-estimates mapping using `retm_estimate`.
3. Re-runs `anc_system.m` with ReTM map.
4. Displays combined comparison plots.
5. Saves comparison data to:
	 - `Simulation_ANC/rm_retm_comparison_results.mat`

## Expected Outputs

- MATLAB figures for time-domain and frequency-domain ANC performance.
- Console summary including average noise reduction in selected bands.
- Optional `.mat` result archive for RM vs ReTM comparison.

## Common Issues

- Path errors:
	- Ensure current folder is `Simulation_ANC` before running scripts.
- Missing MEX/function errors for RIR:
	- Confirm `RIR_Generator/rir_generator.mexw64` is present and loadable.
- Missing source audio:
	- Verify required `.wav` files exist under `Simulation_ANC/Sound_sources/`.

## Suggested Workflow

1. Run `new_anc_simulation.m` first to validate baseline ANC flow.
2. Run `compare_rm_retm_psd.m` to evaluate RM vs ReTM differences.
3. Inspect `rm_retm_comparison_results.mat` for post-analysis or reporting.
