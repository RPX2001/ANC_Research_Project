# RM (Room Model) Estimation Documentation

## Table of Contents
1. [Overview](#overview)
2. [Mathematical Background](#mathematical-background)
3. [Implementation Details](#implementation-details)
4. [Usage Guide](#usage-guide)
5. [Comparison with ReTM and ReTF](#comparison-with-retm-and-retf)
6. [Code Structure](#code-structure)
7. [Examples](#examples)

---

## Overview

### What is RM Estimation?

**RM (Room Model)** estimation is a signal processing technique used to estimate the frequency-domain transfer function mapping between multi-channel reference signals and multi-channel target signals. The RM represents how acoustic signals propagate from reference microphone positions to target microphone positions in a room environment.

**Purpose**: RM estimation is commonly used in:
- Active Noise Control (ANC) systems
- Acoustic echo cancellation
- Spatial audio processing
- Multi-channel noise reduction
- Room acoustic analysis

### Key Characteristics

- **Multi-channel to multi-channel mapping**: Maps B reference channels to A target channels
- **Frequency-domain approach**: Operates in the frequency domain using STFT
- **Least squares estimation**: Uses regularized least squares for robust estimation
- **Assumes linear time-invariant (LTI) system**: The room acoustics are assumed to be linear and time-invariant during the estimation period

---

## Mathematical Background

### Problem Formulation

Given:
- **Target signals**: $\mathbf{e}(t) \in \mathbb{R}^{A}$ (A channels)
- **Reference signals**: $\mathbf{m}(t) \in \mathbb{R}^{B}$ (B channels)

We want to estimate the Room Model **RM** such that:

$$\mathbf{E}(f) = \mathbf{RM}(f) \cdot \mathbf{M}(f)$$

Where:
- $\mathbf{E}(f)$ is the STFT of target signals: $[F \times T \times A]$
- $\mathbf{M}(f)$ is the STFT of reference signals: $[F \times T \times B]$
- $\mathbf{RM}(f)$ is the frequency-domain room model: $[F \times 1 \times A \times B]$

### Least Squares Solution

For each frequency bin $f$, we solve a multi-channel ridge regression problem:

$$\min_{\mathbf{W}_f} \|\mathbf{E}_f - \mathbf{M}_f \mathbf{W}_f^T\|_F^2$$

Where:
- $\mathbf{M}_f \in \mathbb{C}^{T \times B}$ (reference STFT at frequency $f$)
- $\mathbf{E}_f \in \mathbb{C}^{T \times A}$ (target STFT at frequency $f$)
- $\mathbf{W}_f \in \mathbb{C}^{B \times A}$ (filter coefficients at frequency $f$)

The regularized least squares solution is:

$$\mathbf{W}_f = (\mathbf{M}_f^H \mathbf{M}_f + \lambda \mathbf{I})^{-1} \mathbf{M}_f^H \mathbf{E}_f$$

Where:
- $\mathbf{M}_f^H$ is the conjugate transpose (Hermitian) of $\mathbf{M}_f$
- $\lambda$ is the regularization parameter (default: $10^{-6}$)
- $\mathbf{I}$ is the identity matrix of size $B \times B$

The estimated RM at frequency $f$ is:

$$\mathbf{RM}(f) = \mathbf{W}_f^T$$

### Regularization

The regularization term $\lambda \mathbf{I}$ is added to:
1. **Prevent matrix singularity**: Ensures the Gram matrix is invertible
2. **Improve numerical stability**: Reduces sensitivity to noise
3. **Control overfitting**: Prevents overfit to training data

---

## Implementation Details

### File Location
```
Simulation_ANC/ReTM/rm_estimate.m
```

### Function Signature

```matlab
function [rm, rmi] = rm_estimate(tgt_sig, ref_sig, options)
```

### Inputs

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tgt_sig` | `[samples, A]` | Required | Target signals (A channels) |
| `ref_sig` | `[samples, B]` | Required | Reference signals (B channels) |
| `options.fs` | scalar | `[]` | Sampling frequency (Hz) |
| `options.wlen` | scalar | `fs/1` or `4096` | STFT window length (samples) |
| `options.nfft` | scalar | `2^nextpow2(2*wlen)` | FFT size (samples) |
| `options.hop` | scalar | `wlen/4` | STFT hop size (samples) |
| `options.reg` | scalar | `1e-6` | Regularization parameter |

### Outputs

| Output | Type | Description |
|--------|------|-------------|
| `rm` | `[f, 1, A, B]` | Frequency-domain RM (complex) |
| `rmi` | `[samples, 1, A, B]` | Time-domain impulse responses (real) |

### Algorithm Steps

#### 1. **Parameter Initialization**
```matlab
% Set default parameters
wlen = fs / 1;  % 1 second window if fs provided, else 4096
nfft = 2^nextpow2(2 * wlen);  % Zero-padding for better frequency resolution
hop = wlen / 4;  % 75% overlap
reg = 1e-6;  % Regularization
```

#### 2. **STFT Computation**
```matlab
% Create COLA window
wind = shaasp.cola_window(wlen, hop, 'ola');

% Compute STFT with single-sided spectrum
E = shaasp.cola_stft(tgt_sig, nfft, hop, wind, true);  % [F, T, A]
M = shaasp.cola_stft(ref_sig, nfft, hop, wind, true);  % [F, T, B]
```

- Uses Constant Overlap-Add (COLA) windowing
- Single-sided spectrum for computational efficiency
- Results in `F = nfft/2 + 1` frequency bins

#### 3. **Frequency-wise Least Squares**
```matlab
for f = 1:Fbins
    Mf = squeeze(M(f,:,:));  % [T x B]
    Ef = squeeze(E(f,:,:));  % [T x A]
    
    % Regularized least squares
    G = (Mf' * Mf) + reg * eye(B);
    W = G \ (Mf' * Ef);  % [B x A]
    
    rm(f,1,:,:) = W.';  % Store transposed [A x B]
end
```

For each frequency bin:
- Extract time-frequency data for that bin
- Form Gram matrix: $\mathbf{G} = \mathbf{M}_f^H \mathbf{M}_f + \lambda \mathbf{I}$
- Solve: $\mathbf{W}_f = \mathbf{G}^{-1} \mathbf{M}_f^H \mathbf{E}_f$
- Store transposed filter coefficients

#### 4. **Inverse FFT for Impulse Response**
```matlab
rmi = ifft(rm, nfft, 1, 'symmetric');
```

- Converts frequency-domain RM to time-domain impulse responses
- Uses `'symmetric'` flag to ensure real-valued output
- Results in room impulse responses (RIRs) for each channel pair

---

## Usage Guide

### Basic Usage

```matlab
% Load or generate signals
tgt_sig = [samples x A_channels];  % Target microphones
ref_sig = [samples x B_channels];  % Reference microphones

% Estimate RM
[rm, rmi] = rm_estimate(tgt_sig, ref_sig);

% rm  -> [frequencies x 1 x A x B] frequency-domain
% rmi -> [samples x 1 x A x B] time-domain impulse responses
```

### With Parameters

```matlab
% Specify STFT parameters
fs = 16000;  % Sampling rate

[rm, rmi] = rm_estimate(tgt_sig, ref_sig, ...
    'fs', fs, ...
    'wlen', fs/2, ...      % 0.5 second window
    'nfft', 2^15, ...      % 32768 points FFT
    'hop', fs/8, ...       % 0.125 second hop
    'reg', 1e-5);          % Regularization
```

### Apply RM for Filtering

```matlab
% Estimate RM from training data
[rm, ~] = rm_estimate(tgt_train, ref_train, 'fs', fs);

% Apply to new reference signals
ref_new_stft = shaasp.cola_stft(ref_new, nfft, hop, wind, true);

% Multiply each frequency bin
est_tgt_stft = zeros(size(ref_new_stft,1), size(ref_new_stft,2), A);
for f = 1:size(rm,1)
    for t = 1:size(ref_new_stft,2)
        rm_f = squeeze(rm(f,1,:,:));      % [A x B]
        ref_f = squeeze(ref_new_stft(f,t,:));  % [B x 1]
        est_tgt_stft(f,t,:) = rm_f * ref_f;    % [A x 1]
    end
end

% Inverse STFT to time domain
est_tgt = shaasp.cola_istft(est_tgt_stft, nfft, hop, wind);
```

### Noise Reduction Example

```matlab
% Use RM to estimate and remove noise
% Assume: tgt = signal + noise, ref = noise_only

% 1. Estimate RM during noise-only period
[rm, ~] = rm_estimate(tgt_noise, ref_noise, 'fs', fs);

% 2. Estimate noise in target from reference
noise_est_stft = apply_rm_filter(rm, ref_signal_stft);

% 3. Subtract noise from target
clean_stft = tgt_signal_stft - noise_est_stft;

% 4. Convert back to time domain
clean_signal = shaasp.cola_istft(clean_stft, nfft, hop, wind);
```

---

## Comparison with ReTM and ReTF

The repository contains three related estimation methods:

### 1. RM (Room Model) - `rm_estimate.m`

**Approach**: Frequency-wise regularized least squares

**Mathematical Form**:
$$\mathbf{W}_f = (\mathbf{M}_f^H \mathbf{M}_f + \lambda \mathbf{I})^{-1} \mathbf{M}_f^H \mathbf{E}_f$$

**Characteristics**:
- ✅ Simple and fast implementation
- ✅ Direct frequency-domain filtering
- ✅ Explicit regularization control
- ⚠️ Uses direct matrix inversion (may be less stable)
- ⚠️ Frame-by-frame least squares (doesn't use spatial correlation)

**Best for**: 
- Quick prototyping
- Direct filtering applications
- When you need explicit control over regularization

---

### 2. ReTM (Relative Transfer Matrix) - `retm_estimate.m`

**Approach**: Spatial correlation matrix estimation

**Mathematical Form**:
$$\mathbf{ReTM}(f) = \mathbf{R}_{AA}(f) \cdot \mathbf{R}_{BA}^{-1}(f)$$

Where:
- $\mathbf{R}_{AA}$ = Auto-correlation of target signals
- $\mathbf{R}_{BA}$ = Cross-correlation between reference and target

**Characteristics**:
- ✅ Uses correlation matrices (more statistically robust)
- ✅ Better for noisy environments
- ✅ Accounts for spatial relationships
- ✅ Uses pseudo-inverse for better numerical stability
- ⚠️ More computationally intensive
- ⚠️ Requires time-averaging of cross-spectra

**Best for**:
- Multi-channel noise reduction
- Spatial filtering
- Reverberant environments
- Production ANC systems

---

### 3. ReTF (Relative Transfer Function) - `retf_estimate.m`

**Approach**: Single-channel correlation ratio

**Mathematical Form**:
$$\text{ReTF}(f) = \frac{\Phi_{12}(f)}{\Phi_{22}(f)}$$

Where:
- $\Phi_{12}$ = Cross-power spectral density between target and reference
- $\Phi_{22}$ = Auto-power spectral density of reference

**Characteristics**:
- ✅ Simplest implementation
- ✅ Very fast computation
- ✅ Good for single-channel scenarios
- ⚠️ Only works for single-channel (1 target, 1 reference)
- ⚠️ No spatial information

**Best for**:
- Single microphone pairs
- Real-time applications
- Simple noise cancellation

---

### Summary Comparison Table

| Feature | RM | ReTM | ReTF |
|---------|-----|------|------|
| **Channels** | Multi → Multi | Multi → Multi | Single → Single |
| **Method** | Least Squares | Correlation Matrix | PSD Ratio |
| **Complexity** | Medium | High | Low |
| **Robustness** | Medium | High | Low |
| **Speed** | Fast | Moderate | Very Fast |
| **Regularization** | Explicit | Implicit (pinv) | None |
| **Use Case** | General filtering | Advanced ANC | Simple scenarios |

---

## Code Structure

### Dependencies

The `rm_estimate.m` function depends on the `+shaasp` package:

```
+shaasp/
├── cola_window.m          - Generate COLA-compliant window
├── cola_stft.m            - Short-Time Fourier Transform
├── cola_istft.m           - Inverse STFT
└── spec_to_crossspec.m    - Cross-spectral density (used by ReTM)
```

### Related Functions

```
ReTM/
├── rm_estimate.m           - This implementation (RM)
├── retm_estimate.m         - Correlation-based multi-channel
├── retf_estimate.m         - PSD-based single-channel
├── example_retm_denoise.m  - Example usage of ReTM
└── example_retf_denoise.m  - Example usage of ReTF
```

---

## Examples

### Example 1: Basic RM Estimation

```matlab
% Generate test signals
fs = 16000;
duration = 5;  % seconds
t = (0:1/fs:duration-1/fs)';

% Simulate 2 reference sources
ref1 = randn(length(t), 1);
ref2 = randn(length(t), 1);
ref_sig = [ref1, ref2];

% Simulate 3 target signals as filtered versions
h1 = [1; 0.5; 0.2];  % Simple FIR filter
h2 = [0.8; 0.4];
tgt1 = filter(h1, 1, ref1) + filter(h2, 1, ref2);
tgt2 = filter([0.6; 0.3], 1, ref1) + filter([0.9; 0.1], 1, ref2);
tgt3 = filter([0.4; 0.2; 0.1], 1, ref1) + filter([0.7], 1, ref2);
tgt_sig = [tgt1, tgt2, tgt3];

% Add some noise
tgt_sig = tgt_sig + 0.01 * randn(size(tgt_sig));

% Estimate RM
[rm, rmi] = rm_estimate(tgt_sig, ref_sig, 'fs', fs);

% Display results
fprintf('RM size (freq-domain): [%s]\n', num2str(size(rm)));
fprintf('RMI size (time-domain): [%s]\n', num2str(size(rmi)));

% Plot impulse response for target 1, reference 1
figure;
plot(squeeze(rmi(:,1,1,1)));
title('Room Impulse Response: Target 1 ← Reference 1');
xlabel('Samples');
ylabel('Amplitude');
```

### Example 2: Noise Cancellation Application

```matlab
% Scenario: Remove fan noise from speech recording
fs = 16000;

% Load signals
% tgt_with_noise = speech + fan_noise at target mic
% ref_fan_noise = fan noise at reference mic

% Training phase: Estimate RM during speech-free period
noise_only_start = 1;
noise_only_end = 3 * fs;  % First 3 seconds
tgt_noise = tgt_with_noise(noise_only_start:noise_only_end, :);
ref_noise = ref_fan_noise(noise_only_start:noise_only_end, :);

[rm, ~] = rm_estimate(tgt_noise, ref_noise, ...
    'fs', fs, ...
    'wlen', fs/2, ...  % 0.5 sec window for stationary noise
    'reg', 1e-5);

% Apply RM to full reference signal
wlen = fs/2;
nfft = 2^nextpow2(2*wlen);
hop = wlen/4;
wind = shaasp.cola_window(wlen, hop, 'ola');

ref_full_stft = shaasp.cola_stft(ref_fan_noise, nfft, hop, wind, true);
tgt_full_stft = shaasp.cola_stft(tgt_with_noise, nfft, hop, wind, true);

% Estimate noise component in target
[F, T, B] = size(ref_full_stft);
A = size(tgt_full_stft, 3);
noise_est_stft = zeros(F, T, A);

for f = 1:F
    rm_f = squeeze(rm(f,1,:,:));  % [A x B]
    for t = 1:T
        ref_ft = squeeze(ref_full_stft(f,t,:));  % [B x 1]
        noise_est_stft(f,t,:) = rm_f * ref_ft;
    end
end

% Subtract noise from target
clean_stft = tgt_full_stft - noise_est_stft;

% Convert back to time domain
clean_signal = shaasp.cola_istft(clean_stft, nfft, hop, wind);

% Save result
audiowrite('denoised_speech.wav', clean_signal, fs);
```

### Example 3: Analyzing Room Acoustics

```matlab
% Estimate and analyze room impulse responses
fs = 48000;

% tgt_sig = microphone array recordings
% ref_sig = loudspeaker signal

[~, rmi] = rm_estimate(tgt_sig, ref_sig, 'fs', fs, 'wlen', fs);

% Extract and analyze RIRs
A = size(tgt_sig, 2);
B = size(ref_sig, 2);

figure;
for a = 1:A
    for b = 1:B
        subplot(A, B, (a-1)*B + b);
        rir = squeeze(rmi(:,1,a,b));
        
        % Plot impulse response
        t = (0:length(rir)-1) / fs * 1000;  % milliseconds
        plot(t, rir);
        title(sprintf('Target %d ← Ref %d', a, b));
        xlabel('Time (ms)');
        ylabel('Amplitude');
        grid on;
        
        % Calculate RT60 (if needed)
        % rt60 = shaasp.basic_rt60_from_impulse(rir, fs);
        % fprintf('RT60 (T%d←R%d): %.2f ms\n', a, b, rt60*1000);
    end
end
```

---

## Best Practices

### 1. **Choose Appropriate Window Length**
- Should be longer than the room's RT60 (reverberation time)
- Typical values: 0.5 - 2 seconds for rooms
- Longer windows → better frequency resolution, worse time resolution

### 2. **Regularization Parameter Selection**
- Default `1e-6` works for most cases
- Increase for poorly conditioned systems or noisy signals
- Decrease for very clean, well-conditioned data

### 3. **Training Data Requirements**
- Use signals where both target and reference contain the same sources
- Longer training duration → more reliable estimates
- Ensure stationarity during training period

### 4. **Computational Efficiency**
- Use `single_sided = true` in STFT to reduce computation
- Pre-allocate output arrays
- Consider using ReTF for single-channel scenarios

### 5. **Validation**
- Check condition number of Gram matrix: `cond(Mf' * Mf)`
- Verify causality of impulse responses
- Cross-validate with held-out data

---

## Advanced Topics

### Time-Varying RM

For non-stationary environments, estimate RM using shorter time windows:

```matlab
% Sliding window RM estimation
window_duration = 10;  % seconds
window_samples = window_duration * fs;
hop_duration = 5;  % seconds
hop_samples = hop_duration * fs;

n_windows = floor((length(tgt_sig) - window_samples) / hop_samples) + 1;
rm_time_varying = cell(n_windows, 1);

for w = 1:n_windows
    start_idx = (w-1) * hop_samples + 1;
    end_idx = start_idx + window_samples - 1;
    
    tgt_window = tgt_sig(start_idx:end_idx, :);
    ref_window = ref_sig(start_idx:end_idx, :);
    
    rm_time_varying{w} = rm_estimate(tgt_window, ref_window, 'fs', fs);
end
```

### Adaptive Regularization

Adjust regularization based on frequency bin:

```matlab
% Lower frequencies may need more regularization
reg_per_freq = logspace(-7, -5, Fbins);  % Varying regularization

% Modify the loop in rm_estimate.m:
% G = (Mf' * Mf) + reg_per_freq(f) * eye(B);
```

---

## References

1. Birnie, L. (2024). ReTM/ReTF for Acoustic Noise Reduction. Australian National University.
2. Elko, G. W., & Pong, A. T. N. (1995). A simple adaptive first-order differential microphone. IEEE ASSP Workshop on Applications of Signal Processing to Audio and Acoustics.
3. Doclo, S., & Moonen, M. (2002). GSVD-based optimal filtering for single and multimicrophone speech enhancement. IEEE Transactions on Signal Processing.

---

## Contact

For questions or issues regarding the RM estimation implementation:
- **Author**: Lachlan Birnie
- **Organization**: Audio & Acoustic Signal Processing Group - Australian National University
- **Email**: Lachlan.Birnie@anu.edu.au
- **GitHub**: https://github.com/lachlanbirnie

---

*Last Updated: February 19, 2026*
