clearvars;
close all;
clc;

%% Compare RM vs ReTM virtual sensing using existing ANC pipeline
% This script reuses the existing algorithms in:
%   1) new_anc_simulation.m  (tuning + RM estimate + ANC run)
%   2) anc_system.m          (ANC control + PSD/NR computation)
% No ANC algorithm is rewritten here; only estimator selection is switched.

fprintf('=== RM vs ReTM Comparison (same scenario) ===\n');

%% Ensure paths are available
addpath(genpath('RIR_Generator'));
addpath(genpath('ReTM'));

%% 1) Run common tuning stage (unchanged)
fprintf('\n[1/4] Running common tuning stage using new_anc_simulation.m ...\n');
run('new_anc_simulation.m');

if ~exist('rm_est','var') || ~exist('err_sig','var') || ~exist('mon_sig','var') || ~exist('fs','var') || ~exist('wlen','var') || ~exist('nfft','var') || ~exist('hop','var')
    error('Required tuning-stage variables are missing after new_anc_simulation.m.');
end

%% 2) Replace only control-stage primary sources with buccaneer signals
fprintf('[2/4] Preparing control-stage primary sources: buccaneer1.wav + buccaneer2.wav ...\n');

src1_file = fullfile('Sound_sources', 'buccaneer1.wav');
src2_file = fullfile('Sound_sources', 'buccaneer2.wav');
if ~isfile(src1_file)
    error('Control-stage source file not found: %s', src1_file);
end
if ~isfile(src2_file)
    error('Control-stage source file not found: %s', src2_file);
end

[x1_ctrl, fs1_ctrl] = audioread(src1_file);
[x2_ctrl, fs2_ctrl] = audioread(src2_file);

if size(x1_ctrl,2) > 1
    x1_ctrl = mean(x1_ctrl, 2);
end
if size(x2_ctrl,2) > 1
    x2_ctrl = mean(x2_ctrl, 2);
end

if fs1_ctrl ~= fs
    x1_ctrl = resample(x1_ctrl, fs, fs1_ctrl);
end
if fs2_ctrl ~= fs
    x2_ctrl = resample(x2_ctrl, fs, fs2_ctrl);
end

if exist('Ns_target','var')
    N_control = Ns_target;
else
    N_control = size(src,1);
end

if length(x1_ctrl) >= N_control
    x1_ctrl = x1_ctrl(1:N_control);
else
    x1_ctrl = repmat(x1_ctrl, ceil(N_control/length(x1_ctrl)), 1);
    x1_ctrl = x1_ctrl(1:N_control);
end

if length(x2_ctrl) >= N_control
    x2_ctrl = x2_ctrl(1:N_control);
else
    x2_ctrl = repmat(x2_ctrl, ceil(N_control/length(x2_ctrl)), 1);
    x2_ctrl = x2_ctrl(1:N_control);
end

src_control = [x1_ctrl, x2_ctrl];

%% 3) RM control-stage ANC run with new source
fprintf('[3/4] Running ANC control (RM) with buccaneer control-stage sources ...\n');

src = src_control;
run('anc_system.m');

if ~exist('fpsd','var') || ~exist('Pb_dB','var') || ~exist('NR','var') || ~exist('eval_off','var') || ~exist('eval_on','var')
    error('Expected ANC outputs not found after RM control-stage run.');
end

resRM = struct();
resRM.fpsd = fpsd;
resRM.Pb_dB = Pb_dB;
if exist('Pa_dB_shifted','var')
    resRM.Pa_dB = Pa_dB_shifted;
elseif exist('Pa_dB','var')
    resRM.Pa_dB = Pa_dB;
else
    error('ANC ON PSD variable not found for RM run.');
end
resRM.NR = NR;
resRM.eval_off = eval_off;
resRM.eval_on = eval_on;
resRM.t = t;
resRM.pPlot = pPlot;

%% 4) ReTM estimate, then rerun same ANC control script with same source
fprintf('[4/4] Switching estimator to ReTM and rerunning anc_system.m ...\n');

if ~exist('err_sig','var') || ~exist('mon_sig','var') || ~exist('fs','var') || ~exist('wlen','var') || ~exist('nfft','var') || ~exist('hop','var')
    error('Required tuning variables for ReTM estimation are missing.');
end

[retm_est, reim_est] = retm_estimate( ...
    err_sig, ...
    mon_sig, ...
    'fs', fs, ...
    'wlen', wlen, ...
    'nfft', nfft, ...
    'hop', hop);

% anc_system.m currently reads variable rm_est for virtual mapping.
rm_est = retm_est;

% Keep same control-stage source for fair RM vs ReTM comparison
src = src_control;

run('anc_system.m');

if ~exist('fpsd','var') || ~exist('Pb_dB','var') || ~exist('NR','var') || ~exist('eval_off','var') || ~exist('eval_on','var')
    error('Expected ANC outputs not found after ReTM run in anc_system.m.');
end

resReTM = struct();
resReTM.fpsd = fpsd;
resReTM.Pb_dB = Pb_dB;
if exist('Pa_dB_shifted','var')
    resReTM.Pa_dB = Pa_dB_shifted;
elseif exist('Pa_dB','var')
    resReTM.Pa_dB = Pa_dB;
else
    error('ANC ON PSD variable not found for ReTM run.');
end
resReTM.NR = NR;
resReTM.eval_off = eval_off;
resReTM.eval_on = eval_on;
resReTM.t = t;
resReTM.pPlot = pPlot;

%% 3) Plot both methods together (same diagrams)
fprintf('[3/3] Plotting combined PSD and noise reduction curves ...\n');

% Align frequency grids if needed
if length(resRM.fpsd) == length(resReTM.fpsd) && all(abs(resRM.fpsd - resReTM.fpsd) < 1e-12)
    f = resRM.fpsd;
    Pb_rm = resRM.Pb_dB;
    Pa_rm = resRM.Pa_dB;
    NR_rm = resRM.NR;
    Pb_retm = resReTM.Pb_dB;
    Pa_retm = resReTM.Pa_dB;
    NR_retm = resReTM.NR;
else
    f = resRM.fpsd;
    Pb_rm = resRM.Pb_dB;
    Pa_rm = resRM.Pa_dB;
    NR_rm = resRM.NR;
    Pb_retm = interp1(resReTM.fpsd, resReTM.Pb_dB, f, 'linear', 'extrap');
    Pa_retm = interp1(resReTM.fpsd, resReTM.Pa_dB, f, 'linear', 'extrap');
    NR_retm = interp1(resReTM.fpsd, resReTM.NR, f, 'linear', NaN);
end

% Combined PSD: evaluation mic ANC OFF/ON for both RM and ReTM
figure('Name', sprintf('RM vs ReTM - Evaluation Mic %d PSD', resRM.pPlot));
plot(f, Pb_rm, 'k-', 'LineWidth', 1.5); hold on;
plot(f, Pa_rm, 'b-', 'LineWidth', 1.5);
plot(f, Pb_retm, 'k--', 'LineWidth', 1.5);
plot(f, Pa_retm, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
legend('RM: ANC OFF', 'RM: ANC ON', 'ReTM: ANC OFF', 'ReTM: ANC ON', 'Location', 'best');
title(sprintf('Evaluation Microphone %d - PSD Comparison', resRM.pPlot));
xlim([0 600]);

% Combined Noise Reduction
figure('Name', sprintf('RM vs ReTM - Evaluation Mic %d Noise Reduction', resRM.pPlot));
plot(f, NR_rm, 'b', 'LineWidth', 1.8); hold on;
plot(f, NR_retm, 'r', 'LineWidth', 1.8);
grid on;
xlabel('Frequency (Hz)');
ylabel('Noise Reduction (dB)');
legend('RM', 'ReTM', 'Location', 'best');
title(sprintf('Evaluation Microphone %d - Noise Reduction Comparison', resRM.pPlot));
xlim([0 600]);

% Optional: time-domain comparison at evaluation mic
figure('Name', sprintf('RM vs ReTM - Evaluation Mic %d Time Domain', resRM.pPlot));
plot(resRM.t, resRM.eval_off(:,resRM.pPlot), 'k-', 'LineWidth', 1.0); hold on;
plot(resRM.t, resRM.eval_on(:,resRM.pPlot), 'b-', 'LineWidth', 1.0);
plot(resReTM.t, resReTM.eval_on(:,resReTM.pPlot), 'r-', 'LineWidth', 1.0);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
legend('ANC OFF (baseline)', 'ANC ON - RM', 'ANC ON - ReTM', 'Location', 'best');
title(sprintf('Evaluation Microphone %d - Time Domain Comparison', resRM.pPlot));
xlim([resRM.t(1) resRM.t(end)]);

% Summary metrics
fband_low = f >= 50 & f <= 600;
fband_wide = f >= 50 & f <= 1500;

NR_rm_low = nanmean(NR_rm(fband_low));
NR_rm_wide = nanmean(NR_rm(fband_wide));
NR_retm_low = nanmean(NR_retm(fband_low));
NR_retm_wide = nanmean(NR_retm(fband_wide));

fprintf('\n=== RM vs ReTM Summary ===\n');
fprintf('RM   avg NR (50-600 Hz):  %.2f dB\n', NR_rm_low);
fprintf('ReTM avg NR (50-600 Hz):  %.2f dB\n', NR_retm_low);
fprintf('RM   avg NR (50-1500 Hz): %.2f dB\n', NR_rm_wide);
fprintf('ReTM avg NR (50-1500 Hz): %.2f dB\n', NR_retm_wide);

% Save comparison outputs for later use
save('rm_retm_comparison_results.mat', 'resRM', 'resReTM', 'f', 'Pb_rm', 'Pa_rm', 'Pb_retm', 'Pa_retm', 'NR_rm', 'NR_retm');

fprintf('Saved results to rm_retm_comparison_results.mat\n');
fprintf('=== Comparison complete ===\n');
