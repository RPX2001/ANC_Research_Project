%% ANC Control Stage using ReTM
% This script should be run AFTER new_anc_simulation.m
% It expects these variables to exist in the workspace:
% - h_pm, h_sm, h_pee, h_see (impulse responses)
% - retm_est (ReTM estimate)
% - fs, wlen, hop, nfft, win (STFT parameters)
% - src (source signals)

%% Verify required variables exist
required_vars = {'h_pm', 'h_sm', 'h_pee', 'h_see', 'retm_est', 'fs', 'wlen', 'hop', 'nfft', 'win'};
missing_vars = {};
for i = 1:length(required_vars)
    if ~exist(required_vars{i}, 'var')
        missing_vars{end+1} = required_vars{i};
    end
end
if ~isempty(missing_vars)
    error('Missing required variables: %s\nPlease run new_anc_simulation.m first.', strjoin(missing_vars, ', '));
end

fprintf('\n=== Starting ANC Control Stage ===\n');

%% Parameters for ANC Control Stage

% Define control time window (portion of the signal where ANC is active)
T_duration = size(src,1) / fs;    % Total duration of source signal
t_start = 20;                     % seconds - start ANC control  
t_end   = 80;                     % seconds - end ANC control
n_start = floor(t_start*fs) + 1;
n_end   = floor(t_end*fs);

% Extract the control segment from full source signal
src_full = src;                   % Keep full source for reference
src = src_full(n_start:n_end, :); % Crop to control window

% Determine system dimensions
Ns = size(src,1);                 % Number of samples in control window
K  = size(src,2);                 % Number of reference channels (primary sources)
M  = size(h_pm,2);                % Number of monitoring mics
Lspk = size(h_sm,3);              % Number of secondary speakers
F = nfft/2 + 1;                   % Number of frequency bins
Neval = size(h_pee,2);            % Number of evaluation mics
W_nums = 128;                     % Number of FIR filter taps (increased from 64 for better HF modeling)

fprintf('System dimensions: K=%d sources, M=%d mon mics, L=%d speakers, Neval=%d eval mics\n', ...
    K, M, Lspk, Neval);

% Extract ReTM in correct format [F x V x M]
RVM = squeeze(retm_est(:,1,:,:)); % Assuming retm_est is [F x nFrames x V x M]
V = size(RVM,2);                  % Number of virtual error mics

fprintf('ReTM extracted: %d virtual mics mapped from %d monitoring mics\n', V, M);


%% Build STFT of source signals and frequency responses of acoustic paths

fprintf('Computing STFTs of source signals...\n');

X = cell(K,1);
for k = 1:K
    X{k} = stft(src(:,k), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided');  % [F x nFrames]
end
nFrames = size(X{1},2);
fprintf('  STFT frames: %d\n', nFrames);

% Compute frequency responses of all acoustic paths
fprintf('Computing frequency responses of acoustic paths...\n');

% Primary to monitoring: Hpm [F x M x K]
Hpm  = zeros(F, M, K);
Hsm  = zeros(F, M, Lspk);
Hpee = zeros(F, Neval, K);
Hsee = zeros(F, Neval, Lspk);

for k = 1:K
    for m = 1:M
        H = fft(h_pm(:,m,k), nfft);
        Hpm(:,m,k) = H(1:F);
    end
    for p = 1:Neval
        H = fft(h_pee(:,p,k), nfft);
        Hpee(:,p,k) = H(1:F);
    end
end

for l = 1:Lspk
    for m = 1:M
        H = fft(h_sm(:,m,l), nfft);
        Hsm(:,m,l) = H(1:F);
    end
    for p = 1:Neval
        H = fft(h_see(:,p,l), nfft);
        Hsee(:,p,l) = H(1:F);
    end
end


%% Compute primary disturbance at monitoring and evaluation mics (baseline - ANC OFF)

fprintf('Computing primary disturbance signals...\n');

% Primary contribution at monitoring mics: Dm_tf [F x nFrames x M]
Dm_tf = zeros(F, nFrames, M);
for k = 1:K
    for m = 1:M
        Dm_tf(:,:,m) = Dm_tf(:,:,m) + Hpm(:,m,k) .* X{k};
    end
end

% Primary contribution at evaluation mics (ANC OFF baseline): Deval_tf_off [F x nFrames x Neval]  
Deval_tf_off = zeros(F, nFrames, Neval);
for k = 1:K
    for p = 1:Neval
        Deval_tf_off(:,:,p) = Deval_tf_off(:,:,p) + Hpee(:,p,k) .* X{k};
    end
end

fprintf('  Primary disturbance STFT computed\n');


%% Initialize ANC adaptive filter and control parameters

fprintf('\n=== Initializing ANC Controller ===\n');

% Adaptive filter W(f,l,k) in frequency domain [F x Lspk x K]
Wf = zeros(F, Lspk, K);           % Frequency-domain filter weights

% Storage for results: virtual error and evaluation mic signals
EV0_tf = zeros(F, nFrames, V);             % Virtual error (ANC OFF) 
EV_tf  = zeros(F, nFrames, V);             % Virtual error (ANC ON)
Deval_tf_on = zeros(F, nFrames, Neval);    % Evaluation mics (ANC ON)

% FxLMS algorithm parameters
mu = 8e-4;                        % Normalized FxLMS step size
pwr_floor = 1e-6;                 % Power floor to prevent division by zero
leak = 1e-5;                      % Small leakage for robustness

% Frequency adaptation range
f_axis = (0:F-1)' * fs/nfft;
adapt_bins = find(f_axis >= 20 & f_axis <= 1500);  % Match paper-like control band
% Optional: limit to specific frequency range if needed
% f_axis = (0:F-1)' * fs/nfft;
% adapt_bins = find(f_axis >= 50 & f_axis <= 1500);  % Focus on specific band

fprintf('  Step size mu = %.1e\n', mu);
fprintf('  Adapting %d frequency bins (%.1f - %.1f Hz)\n', ...
    length(adapt_bins), min(adapt_bins)*fs/nfft, max(adapt_bins)*fs/nfft);


%% Main ANC Control Loop (FxLMS with ReTM virtual sensing)

fprintf('\n=== Running ANC Control Loop ===\n');

% Prepare reference signal matrix: Xmat [F x nFrames x K]
Xmat = zeros(F, nFrames, K);
for k = 1:K
    Xmat(:,:,k) = X{k};
end
K = size(Xmat,3);  % lock to actual reference-channel dimension

% Skip initial and final frames to avoid transient artifacts
skipFrames = 60;
tfrm_start = 1 + skipFrames;
tfrm_end   = nFrames - skipFrames;

fprintf('  Processing frames %d to %d (skipping %d at start/end)\n', ...
    tfrm_start, tfrm_end, skipFrames);

% Progress reporting
report_interval = max(1, floor((tfrm_end - tfrm_start) / 20)); % Report every 5%
tic;

for tfrm = tfrm_start:tfrm_end
    
    % Progress reporting
    if mod(tfrm - tfrm_start, report_interval) == 0
        pct = 100 * (tfrm - tfrm_start) / (tfrm_end - tfrm_start);
        fprintf('  Progress: %.0f%% (frame %d/%d)\n', pct, tfrm, tfrm_end);
    end

    % Reference at this frame: Xf(k) for each bin
    % We compute speaker outputs for ALL frequency bins:
    Yf = zeros(F, Lspk);   % [F x L]

    for l = 1:Lspk
        acc = zeros(F,1);
        for k = 1:K
            acc = acc + Wf(:,l,k) .* Xmat(:,tfrm,k);  % elementwise over F
        end
        Yf(:,l) = acc;
    end

    % Monitoring residual at ALL bins: eM = dM + S_M * y
    eM_tf_frame = squeeze(Dm_tf(:,tfrm,:));   % [F x M]
    for l = 1:Lspk
        eM_tf_frame = eM_tf_frame + squeeze(Hsm(:,:,l)) .* Yf(:,l); % [F x M]
    end

    % Virtual error estimate: eV = RVM * eM  (per frequency bin)
    for fbin = 1:F
        Rf = squeeze(RVM(fbin,:,:));                 % [V x M]

        dMf = squeeze(Dm_tf(fbin,tfrm,:));           % [M x 1]
        EV0_tf(fbin,tfrm,:) = Rf * dMf;              % baseline (ANC OFF)

        eMf = eM_tf_frame(fbin,:).';                 % [M x 1]
        eVf = Rf * eMf;                              % [V x 1]
        EV_tf(fbin,tfrm,:) = eVf;
    end

    % Evaluation mics (ANC ON): Deval = P_eval*X + S_eval*Y
    Deval_frame = squeeze(Deval_tf_off(:,tfrm,:));   % [F x Neval]
    for l = 1:Lspk
        Deval_frame = Deval_frame + squeeze(Hsee(:,:,l)) .* Yf(:,l); % [F x Neval]
    end
    Deval_tf_on(:,tfrm,:) = Deval_frame;

    % --- STFT-domain FxLMS update ---
    for fbin = adapt_bins

        % reference vector
        Xvec = squeeze(Xmat(fbin,tfrm,:));                  % [K x 1] (runtime-safe)
        Xvec = Xvec(:);

        % get error vector at virtual mics
        eVf = squeeze(EV_tf(fbin,tfrm,:));                  % [V x 1] (runtime-safe)
        eVf = eVf(:);

        % Secondary path to monitoring mics (M x L)
        SMf = squeeze(Hsm(fbin,:,:));                       % [M x L]
        SMf = reshape(SMf, size(SMf,1), []);                % force 2-D
        if size(SMf,2) ~= Lspk && size(SMf,1) == Lspk
            SMf = SMf.';                                    % ensure [M x L]
        end

        % ReTM mapping (V x M)
        Rf = squeeze(RVM(fbin,:,:));                        % [V x M]
        Rf = reshape(Rf, size(Rf,1), []);                   % force 2-D
        if size(Rf,1) ~= numel(eVf) && size(Rf,2) == numel(eVf)
            Rf = Rf.';                                      % ensure [V x M]
        end

        % gradient term per speaker:
        % g = S_M^H * R^H * eV
        eVf = reshape(eVf, [], 1);                          % force column vector
        RH  = permute(conj(Rf),  [2 1]);                    % R^H, size [M x V]
        SMH = permute(conj(SMf), [2 1]);                    % S_M^H, size [L x M]

        % Runtime-safe multiply with aligned V-dimension
        V_use = min(size(RH,2), numel(eVf));
        tmp = RH(:,1:V_use) * eVf(1:V_use);                % [M x 1]

        % Runtime-safe multiply in case M dimensions differ by orientation drift
        M_use = min(size(SMH,2), numel(tmp));
        g = SMH(:,1:M_use) * tmp(1:M_use);                 % [L x 1]

        % Runtime-safe channel count (protects against dimension drift)
        K_use = min(size(Wf,3), numel(Xvec));

        % Per-bin normalized step (NLMS style)
        xpow = sum(abs(Xvec(1:K_use)).^2);
        eta = mu / (xpow + pwr_floor);

        % Update weights
        for l = 1:Lspk
            for k = 1:K_use
                Wf(fbin,l,k) = (1 - leak) * Wf(fbin,l,k) - eta * conj(Xvec(k)) * g(l);
            end
        end
    end

end

elapsed = toc;
fprintf('  ANC loop completed in %.1f seconds\n', elapsed);
fprintf('  Average processing speed: %.2f frames/sec\n', (tfrm_end-tfrm_start)/elapsed);

% Fill skipped frames with baseline (ANC OFF) to avoid ISTFT artifacts
fprintf('  Filling skipped frames with baseline signal...\n');
EV_tf(:,1:tfrm_start-1,:) = EV0_tf(:,1:tfrm_start-1,:);
EV_tf(:,tfrm_end+1:end,:) = EV0_tf(:,tfrm_end+1:end,:);
Deval_tf_on(:,1:tfrm_start-1,:) = Deval_tf_off(:,1:tfrm_start-1,:);
Deval_tf_on(:,tfrm_end+1:end,:) = Deval_tf_off(:,tfrm_end+1:end,:);

%% Convert STFT back to time domain (ISTFT)

fprintf('\n=== Converting back to time domain ===\n');

% Helper function to force exact length
force_len = @(x,N) (length(x)>=N)*x(1:N) + (length(x)<N)*[x(:); zeros(N-length(x),1)];

% ISTFT virtual error signals
fprintf('  Converting virtual error signals...\n');
ev_before = zeros(Ns, V);
ev_after  = zeros(Ns, V);
for v = 1:V
    tmp = istft(EV0_tf(:,:,v), fs, 'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided');
    ev_before(:,v) = force_len(tmp, Ns);

    tmp = istft(EV_tf(:,:,v), fs, 'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided');
    ev_after(:,v) = force_len(tmp, Ns);
end

% ISTFT evaluation microphone signals  
fprintf('  Converting evaluation microphone signals...\n');
eval_off = zeros(Ns, Neval);
eval_on  = zeros(Ns, Neval);
for p = 1:Neval
    tmp = istft(Deval_tf_off(:,:,p), fs, 'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided');
    eval_off(:,p) = force_len(tmp, Ns);

    tmp = istft(Deval_tf_on(:,:,p), fs, 'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided');
    eval_on(:,p) = force_len(tmp, Ns);
end

fprintf('  Time-domain signals ready\n');


%% Results Visualization and Analysis

fprintf('\n=== Generating Results ===\n');

% Time axis (relative to full signal: 30-70s corresponds to control window)
t = ((n_start-1) + (0:Ns-1))' / fs;   

% Select which evaluation mic to plot (1 or 2)
pPlot = 1;

%% Plot 1: Time-domain comparison
fprintf('  Creating time-domain plot...\n');
figure('Name', sprintf('Evaluation Mic %d - Time Domain', pPlot));
plot(t, eval_off(:,pPlot), 'b', 'LineWidth', 1); hold on;
plot(t, eval_on(:,pPlot), 'r', 'LineWidth', 1);
grid on; 
xlabel('Time (s)'); 
ylabel('Amplitude');
legend('ANC OFF (baseline)', 'ANC ON', 'Location', 'best');
title(sprintf('Evaluation Microphone %d: Time Domain Comparison', pPlot));
xlim([t(1) t(end)]);  % Show full control window
% Or zoom to specific region: xlim([35 65]);

%% Plot 2: Power Spectral Density (PSD)
fprintf('  Computing PSD...\n');

% PSD computation parameters
nwel = 4096;                      % Welch window length
nover = round(0.75*nwel);         % 75% overlap
nfftW = 4096;                     % FFT length for Welch

% Trim transients at beginning and end
trim = 2*wlen;
xb = eval_off(trim:end-trim, pPlot);  % ANC OFF (baseline)
xa = eval_on(trim:end-trim, pPlot);   % ANC ON

% Optional: use virtual error instead of evaluation mics
% xb = ev_before(trim:end-trim, pPlot);
% xa = ev_after(trim:end-trim, pPlot);


% Compute PSD using Welch method
[Pb, fpsd] = pwelch(xb, hann(nwel,'periodic'), nover, nfftW, fs);
[Pa, ~]    = pwelch(xa, hann(nwel,'periodic'), nover, nfftW, fs);

% Convert to dB scale
Pb_dB = 10*log10(max(Pb, 1e-20));
Pa_dB = 10*log10(max(Pa, 1e-20));

% Compute noise reduction
NR = Pb_dB - Pa_dB;

% Mask low-energy bins to avoid spurious noise reduction values
thr = max(Pb_dB) - 40;            % Threshold: 40 dB below peak
mask = Pb_dB > thr;
NR(~mask) = NaN;                  % Set masked bins to NaN (won't plot)
NR = min(max(NR, -30), 60);       % Clip to reasonable range [-30, 60] dB

fprintf('  PSD computed successfully\n');

% Plot PSD comparison
figure('Name', sprintf('Evaluation Mic %d - PSD', pPlot));
plot(fpsd, Pb_dB, 'b', 'LineWidth', 1.5); hold on; 
plot(fpsd, Pa_dB, 'r', 'LineWidth', 1.5);
grid on; 
xlim([0 fs/2]);
xlabel('Frequency (Hz)'); 
ylabel('PSD (dB/Hz)');
legend('ANC OFF (baseline)', 'ANC ON', 'Location', 'best');
title(sprintf('Evaluation Microphone %d: Power Spectral Density', pPlot));
% Show extended frequency range to visualize high-frequency performance
xlim([0 2000]);

%% Plot 3: Noise Reduction vs Frequency
fprintf('  Creating noise reduction plot...\n');

figure('Name', sprintf('Evaluation Mic %d - Noise Reduction', pPlot));
plot(fpsd, NR, 'LineWidth', 1.5, 'Color', [0 0.5 0]);
grid on; 
xlim([0 fs/2]);
xlabel('Frequency (Hz)'); 
ylabel('Noise Reduction (dB)');
title(sprintf('Evaluation Microphone %d: ANC Performance', pPlot));
xlim([0 2000]);

% Compute average noise reduction in frequency band of interest
fband_low = fpsd >= 50 & fpsd <= 600;     % Traditional ANC band
fband_wide = fpsd >= 50 & fpsd <= 1500;   % Extended band for evaluation
NR_avg_low = nanmean(NR(fband_low));      % Average NR in low band
NR_avg_wide = nanmean(NR(fband_wide));    % Average NR in wide band

fprintf('\n=== Performance Summary ===\n');
fprintf('  Average noise reduction (50-600 Hz): %.1f dB\n', NR_avg_low);
fprintf('  Average noise reduction (50-1500 Hz): %.1f dB\n', NR_avg_wide);
fprintf('  Peak noise reduction: %.1f dB at %.1f Hz\n', ...
    max(NR(fband_wide)), fpsd(find(NR == max(NR(fband_wide)), 1)));

fprintf('\n=== ANC Analysis Complete ===\n');