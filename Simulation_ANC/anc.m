clearvars;
close all;
clc;

fprintf('=== Fig. 6(b) reproduction: Scenario 3 (secondary-path change) ===\n');

addpath(genpath('RIR_Generator'));
addpath(genpath('ReTM'));

rng(0);

%% ============================================================
% Parameters
%% ============================================================
c = 340;
fs = 8000;
Lroom = [8 7 4];
beta = 0.4;
n = 4096;
mtype = 'omnidirectional';
order = -1;
dim = 3;
orientation = 0;
hp_filter = 1;

J = 2;                   % primary sources
K = 2;                   % secondary loudspeakers
M = 8;                   % monitoring microphones
V = 6;                   % virtual microphones
P_eval = 2;              % last pair for evaluation

Tseg = 60;
Ns_seg = fs * Tseg;
N_total = 3 * Ns_seg;

wlen = 1024;
hop  = 256;
nfft = 2048;
win  = hann(wlen,'periodic');
F    = nfft/2 + 1;

%% ============================================================
% Geometry
%% ============================================================
mon_mics = [ 0.15   0.15  -0.15;
             0.15  -0.15  -0.15;
            -0.15   0.15  -0.15;
            -0.15  -0.15  -0.15;
             0      0.212  0.15;
             0     -0.212  0.15;
             0.212  0      0.15;
            -0.212  0      0.15];

prim_sources = [0.6  0.8  1.00;
                2.17 1.15 0.05];

virt_mics = [ 0.05  0.10  0;
              0.05 -0.10  0;
             -0.05  0.10  0;
             -0.05 -0.10  0;
              0     0.10  0;
              0    -0.10  0];

eval_idx = [5 6];
eval_mics = virt_mics(eval_idx,:);

% Baseline tuning secondary positions
sec_sources_tune = [0   0.30 0;
                    0  -0.30 0];

% Scenario 3 changed secondary positions
sec_sources_ctrl = [0   0.35 0;
                    0  -0.35 0];

%% ============================================================
% Load source signals
%% ============================================================
fprintf('Loading source signals...\n');

[x1, fs1] = audioread('Sound_sources/buccaneer1.wav');
[x2, fs2] = audioread('Sound_sources/buccaneer2.wav');

if size(x1,2) > 1, x1 = mean(x1,2); end
if size(x2,2) > 1, x2 = mean(x2,2); end

if fs1 ~= fs, x1 = resample(x1, fs, fs1); end
if fs2 ~= fs, x2 = resample(x2, fs, fs2); end

x1 = make_exact_length(x1, N_total);
x2 = make_exact_length(x2, N_total);

src_all = [x1, x2];

src_tune = src_all(1:Ns_seg, :);
src_scn2 = src_all(Ns_seg+1:2*Ns_seg, :);
src_scn3 = src_all(2*Ns_seg+1:3*Ns_seg, :);

%% ============================================================
% Generate impulse responses
%% ============================================================
fprintf('Generating impulse responses...\n');

h_pm = zeros(n, M, J);
h_pv = zeros(n, V, J);
h_peval = zeros(n, P_eval, J);

for j = 1:J
    h_pm(:,:,j) = rir_generator(c, fs, mon_mics, prim_sources(j,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
    h_pv(:,:,j) = rir_generator(c, fs, virt_mics, prim_sources(j,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
    h_peval(:,:,j) = rir_generator(c, fs, eval_mics, prim_sources(j,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
end

h_sm_tune = zeros(n, M, K);
h_sv_tune = zeros(n, V, K);
h_seval_tune = zeros(n, P_eval, K);

for k = 1:K
    h_sm_tune(:,:,k) = rir_generator(c, fs, mon_mics, sec_sources_tune(k,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
    h_sv_tune(:,:,k) = rir_generator(c, fs, virt_mics, sec_sources_tune(k,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
    h_seval_tune(:,:,k) = rir_generator(c, fs, eval_mics, sec_sources_tune(k,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
end

h_sm_ctrl = zeros(n, M, K);
h_sv_ctrl = zeros(n, V, K);
h_seval_ctrl = zeros(n, P_eval, K);

for k = 1:K
    h_sm_ctrl(:,:,k) = rir_generator(c, fs, mon_mics, sec_sources_ctrl(k,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
    h_sv_ctrl(:,:,k) = rir_generator(c, fs, virt_mics, sec_sources_ctrl(k,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
    h_seval_ctrl(:,:,k) = rir_generator(c, fs, eval_mics, sec_sources_ctrl(k,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
end

%% ============================================================
% Frequency-domain paths
%% ============================================================
H_pm         = rir_to_H(h_pm, nfft);
H_pv         = rir_to_H(h_pv, nfft);
H_peval      = rir_to_H(h_peval, nfft);

H_sm_tune    = rir_to_H(h_sm_tune, nfft);
H_sv_tune    = rir_to_H(h_sv_tune, nfft);
H_seval_tune = rir_to_H(h_seval_tune, nfft);

H_sm_ctrl    = rir_to_H(h_sm_ctrl, nfft);
H_sv_ctrl    = rir_to_H(h_sv_ctrl, nfft);
H_seval_ctrl = rir_to_H(h_seval_ctrl, nfft);

%% ============================================================
% Tuning-stage data
% ReTM: primary + secondary active
% RM  : primary only
%% ============================================================
fprintf('Generating tuning-stage signals for ReTM and RM separately...\n');

X_tune = cell(J,1);
for j = 1:J
    X_tune{j} = stft(src_tune(:,j), fs, ...
        'Window', win, ...
        'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, ...
        'FrequencyRange', 'onesided');
end
nFrames_tune = size(X_tune{1},2);

U_tune = cell(K,1);
sec_noise = randn(Ns_seg, K);
for k = 1:K
    sec_noise(:,k) = bandpass(sec_noise(:,k), [20 600], fs);
    U_tune{k} = stft(sec_noise(:,k), fs, ...
        'Window', win, ...
        'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, ...
        'FrequencyRange', 'onesided');
end

EM_tf_rm = zeros(F, nFrames_tune, M);
EV_tf_rm = zeros(F, nFrames_tune, V);

for j = 1:J
    for m = 1:M
        EM_tf_rm(:,:,m) = EM_tf_rm(:,:,m) + H_pm(:,m,j) .* X_tune{j};
    end
    for v = 1:V
        EV_tf_rm(:,:,v) = EV_tf_rm(:,:,v) + H_pv(:,v,j) .* X_tune{j};
    end
end

EM_tf_retm = EM_tf_rm;
EV_tf_retm = EV_tf_rm;

for k = 1:K
    for m = 1:M
        EM_tf_retm(:,:,m) = EM_tf_retm(:,:,m) + H_sm_tune(:,m,k) .* U_tune{k};
    end
    for v = 1:V
        EV_tf_retm(:,:,v) = EV_tf_retm(:,:,v) + H_sv_tune(:,v,k) .* U_tune{k};
    end
end

mon_sig_rm    = zeros(Ns_seg, M);
virt_sig_rm   = zeros(Ns_seg, V);
mon_sig_retm  = zeros(Ns_seg, M);
virt_sig_retm = zeros(Ns_seg, V);

for m = 1:M
    mon_sig_rm(:,m) = force_len(istft(EM_tf_rm(:,:,m), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided'), Ns_seg);

    mon_sig_retm(:,m) = force_len(istft(EM_tf_retm(:,:,m), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided'), Ns_seg);
end

for v = 1:V
    virt_sig_rm(:,v) = force_len(istft(EV_tf_rm(:,:,v), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided'), Ns_seg);

    virt_sig_retm(:,v) = force_len(istft(EV_tf_retm(:,:,v), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided'), Ns_seg);
end

%% ============================================================
% Estimate ReTM and RM
%% ============================================================
fprintf('Estimating ReTM...\n');
[retm_est, ~] = retm_estimate(virt_sig_retm, mon_sig_retm, ...
    'fs', fs, 'wlen', wlen, 'nfft', nfft, 'hop', hop);

fprintf('Estimating RM...\n');
[rm_est, ~] = rm_estimate(virt_sig_rm, mon_sig_rm, ...
    'fs', fs, 'wlen', wlen, 'nfft', nfft, 'hop', hop, 'reg', 1e-6);

fprintf('size(retm_est) = %s\n', mat2str(size(retm_est)));
fprintf('size(rm_est)   = %s\n', mat2str(size(rm_est)));

%% ============================================================
% Control config
%% ============================================================
cfg = struct();
cfg.mu = 8e-4;
cfg.pwr_floor = 1e-6;
cfg.leak = 1e-5;
cfg.f_lo = 20;
cfg.f_hi = 600;
cfg.skipFrames = 60;

%% ============================================================
% Run scenario 2
%% ============================================================
fprintf('\nRunning Scenario 2 (no sec-path change) with ReTM...\n');
out2_ReTM = run_control( ...
    'ReTM', retm_est, ...
    H_pm, H_pv, H_peval, ...
    H_sm_tune, H_sv_tune, H_seval_tune, ...
    H_sm_tune, H_sv_tune, H_seval_tune, ...
    src_scn2, fs, wlen, hop, nfft, win, cfg);

fprintf('\nRunning Scenario 2 (no sec-path change) with RM...\n');
out2_RM = run_control( ...
    'RM', rm_est, ...
    H_pm, H_pv, H_peval, ...
    H_sm_tune, H_sv_tune, H_seval_tune, ...
    H_sm_tune, H_sv_tune, H_seval_tune, ...
    src_scn2, fs, wlen, hop, nfft, win, cfg);

%% ============================================================
% Run scenario 3
% IMPORTANT FIX:
%   ReTM uses UPDATED monitoring secondary path model H_sm_ctrl
%   RM keeps old virtual/model paths
%% ============================================================
fprintf('\nRunning Scenario 3 (changed sec-path) with ReTM...\n');
out3_ReTM = run_control( ...
    'ReTM', retm_est, ...
    H_pm, H_pv, H_peval, ...
    H_sm_ctrl, H_sv_tune, H_seval_tune, ...
    H_sm_ctrl, H_sv_ctrl, H_seval_ctrl, ...
    src_scn3, fs, wlen, hop, nfft, win, cfg);

fprintf('\nRunning Scenario 3 (changed sec-path) with RM...\n');
out3_RM = run_control( ...
    'RM', rm_est, ...
    H_pm, H_pv, H_peval, ...
    H_sm_tune, H_sv_tune, H_seval_tune, ...
    H_sm_ctrl, H_sv_ctrl, H_seval_ctrl, ...
    src_scn3, fs, wlen, hop, nfft, win, cfg);

%% ============================================================
% Plot spectra
%% ============================================================
trim = 2*wlen;
p = 1;

% Scenario 2
x2_off  = mean(out2_ReTM.eval_off(trim:end-trim, :), 2);
x2_retm = mean(out2_ReTM.eval_on(trim:end-trim, :), 2);
x2_rm   = mean(out2_RM.eval_on(trim:end-trim, :), 2);

[P2b, f_plot] = stft_avg_power(x2_off, fs, win, hop, nfft);
[P2r, ~]      = stft_avg_power(x2_retm, fs, win, hop, nfft);
[P2m, ~]      = stft_avg_power(x2_rm, fs, win, hop, nfft);

P2b_dB = 10*log10(max(P2b,1e-20));
P2r_dB = 10*log10(max(P2r,1e-20));
P2m_dB = 10*log10(max(P2m,1e-20));

figure('Name','Fig 6(a) reproduction');
plot(f_plot, P2r_dB, 'k', 'LineWidth', 1.8); hold on;
plot(f_plot, P2m_dB, '--', 'Color', [1 0 1], 'LineWidth', 1.8);
plot(f_plot, P2b_dB, ':b', 'LineWidth', 2.0);
grid on;
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('ReTM', 'RM', 'No ANC', 'Location', 'best');
title('Scenario 2: Without secondary paths change');
xlim([0 600]);

% Scenario 3
x3_off  = mean(out3_ReTM.eval_off(trim:end-trim, :), 2);
x3_retm = mean(out3_ReTM.eval_on(trim:end-trim, :), 2);
x3_rm   = mean(out3_RM.eval_on(trim:end-trim, :), 2);

[P3b, f_plot] = stft_avg_power(x3_off, fs, win, hop, nfft);
[P3r, ~]      = stft_avg_power(x3_retm, fs, win, hop, nfft);
[P3m, ~]      = stft_avg_power(x3_rm, fs, win, hop, nfft);

P3b_dB = 10*log10(max(P3b,1e-20));
P3r_dB = 10*log10(max(P3r,1e-20));
P3m_dB = 10*log10(max(P3m,1e-20));

P3b_dB = smoothdata(P3b_dB, 'movmean', 7);
P3r_dB = smoothdata(P3r_dB, 'movmean', 7);
P3m_dB = smoothdata(P3m_dB, 'movmean', 7);

figure('Name','Fig 6(b) reproduction');
plot(f_plot, P3r_dB, 'k', 'LineWidth', 1.8); hold on;
plot(f_plot, P3m_dB, '--', 'Color', [1 0 1], 'LineWidth', 1.8);
plot(f_plot, P3b_dB, ':b', 'LineWidth', 2.0);
grid on;
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('ReTM', 'RM', 'No ANC', 'Location', 'best');
title('Scenario 3: With secondary paths change');
xlim([0 600]);

% Scenario 3 NR
NR3_retm = P3b_dB - P3r_dB;
NR3_rm   = P3b_dB - P3m_dB;

figure('Name','Scenario 3 NR');
plot(f_plot, NR3_retm, 'k', 'LineWidth', 1.8); hold on;
plot(f_plot, NR3_rm, '--', 'Color', [1 0 1], 'LineWidth', 1.8);
yline(0, 'k:');
grid on;
xlabel('Frequency (Hz)');
ylabel('Noise Reduction (dB)');
legend('ReTM','RM','Location','best');
title('Scenario 3: Noise reduction');
xlim([0 600]);

idx_low = f_plot >= 20 & f_plot < 400;
idx_hi  = f_plot >= 400 & f_plot <= 600;

fprintf('\nScenario 3 average residual difference (RM - ReTM):\n');
fprintf('20-400 Hz  : %.4f dB\n', mean(P3m_dB(idx_low) - P3r_dB(idx_low)));
fprintf('400-600 Hz : %.4f dB\n', mean(P3m_dB(idx_hi)  - P3r_dB(idx_hi)));

%% ============================================================
% Local functions
%% ============================================================

function y = make_exact_length(x, N)
    x = x(:);
    if numel(x) > N
        y = x(1:N);
    else
        y = repmat(x, ceil(N/numel(x)), 1);
        y = y(1:N);
    end
end

function y = force_len(x, N)
    x = x(:);
    if numel(x) >= N
        y = x(1:N);
    else
        y = [x; zeros(N-numel(x),1)];
    end
end

function H = rir_to_H(h, nfft)
    F = nfft/2 + 1;
    [~, nch, nsrc] = size(h);
    H = zeros(F, nch, nsrc);
    for s = 1:nsrc
        for c = 1:nch
            tmp = fft(h(:,c,s), nfft);
            H(:,c,s) = tmp(1:F);
        end
    end
end

function RVM = extract_map(map_est)
    sz = size(map_est);

    if ndims(map_est) == 4
        if sz(2) == 1
            RVM = squeeze(map_est(:,1,:,:));   % [F x V x M]
        else
            error('Unexpected map_est size %s. Adjust extract_map().', mat2str(sz));
        end
    elseif ndims(map_est) == 3
        RVM = map_est;
    else
        error('Unsupported map_est size %s', mat2str(sz));
    end
end

function [Pavg, f] = stft_avg_power(x, fs, win, hop, nfft)
    S = stft(x, fs, ...
        'Window', win, ...
        'OverlapLength', length(win)-hop, ...
        'FFTLength', nfft, ...
        'FrequencyRange', 'onesided');
    P = abs(S).^2;
    Pavg = mean(P, 2);
    f = (0:nfft/2)' * fs/nfft;
end

function out = run_control(methodName, map_est, ...
    H_pm, H_pv, H_peval, ...
    H_sm_model, H_sv_model, H_seval_model, ...
    H_sm_actual, H_sv_actual, H_seval_actual, ...
    src_ctrl, fs, wlen, hop, nfft, win, cfg)

    F = nfft/2 + 1;
    J = size(H_pm, 3);
    M = size(H_pm, 2);
    V = size(H_pv, 2);
    K = size(H_sm_actual, 3);
    P_eval = size(H_peval, 2);
    Ns = size(src_ctrl, 1);

    RVM = extract_map(map_est);   % [F x V x M]

    X = cell(J,1);
    for j = 1:J
        X{j} = stft(src_ctrl(:,j), fs, ...
            'Window', win, ...
            'OverlapLength', wlen-hop, ...
            'FFTLength', nfft, ...
            'FrequencyRange', 'onesided');
    end
    nFrames = size(X{1},2);

    Xmat = zeros(F, nFrames, J);
    for j = 1:J
        Xmat(:,:,j) = X{j};
    end

    Dm_tf     = zeros(F, nFrames, M);
    Dv_tf_off = zeros(F, nFrames, V);
    De_tf_off = zeros(F, nFrames, P_eval);

    for j = 1:J
        for m = 1:M
            Dm_tf(:,:,m) = Dm_tf(:,:,m) + H_pm(:,m,j) .* Xmat(:,:,j);
        end
        for v = 1:V
            Dv_tf_off(:,:,v) = Dv_tf_off(:,:,v) + H_pv(:,v,j) .* Xmat(:,:,j);
        end
        for p = 1:P_eval
            De_tf_off(:,:,p) = De_tf_off(:,:,p) + H_peval(:,p,j) .* Xmat(:,:,j);
        end
    end

    Wf = zeros(F, K, J);
    Dv_tf_on = zeros(F, nFrames, V);
    De_tf_on = zeros(F, nFrames, P_eval);

    f_axis = (0:F-1)' * fs/nfft;
    adapt_bins = find(f_axis >= cfg.f_lo & f_axis <= cfg.f_hi);

    skipFrames = cfg.skipFrames;
    tfrm_start = 1 + skipFrames;
    tfrm_end   = nFrames - skipFrames;

    for tfrm = tfrm_start:tfrm_end

        Yf = zeros(F, K);
        for k = 1:K
            acc = zeros(F,1);
            for j = 1:J
                acc = acc + Wf(:,k,j) .* Xmat(:,tfrm,j);
            end
            Yf(:,k) = acc;
        end

        eM_frame = squeeze(Dm_tf(:,tfrm,:));   % [F x M]
        for k = 1:K
            eM_frame = eM_frame + squeeze(H_sm_actual(:,:,k)) .* Yf(:,k);
        end

        eV_actual = squeeze(Dv_tf_off(:,tfrm,:));  % [F x V]
        for k = 1:K
            eV_actual = eV_actual + squeeze(H_sv_actual(:,:,k)) .* Yf(:,k);
        end
        Dv_tf_on(:,tfrm,:) = eV_actual;

        eEval_actual = squeeze(De_tf_off(:,tfrm,:)); % [F x P_eval]
        for k = 1:K
            eEval_actual = eEval_actual + squeeze(H_seval_actual(:,:,k)) .* Yf(:,k);
        end
        De_tf_on(:,tfrm,:) = eEval_actual;

        for ii = 1:numel(adapt_bins)
            fbin = adapt_bins(ii);

            Xvec = reshape(Xmat(fbin,tfrm,:), [J,1]);
            eMf  = reshape(eM_frame(fbin,:), [M,1]);

            if strcmpi(methodName, 'ReTM')
                Rf = reshape(RVM(fbin,:,:), [V,M]);
                eV_hat = Rf * eMf;

                % IMPORTANT: now uses H_sm_model passed from caller
                SMm = reshape(H_sm_model(fbin,:,:), [M,K]);
                Gf = Rf * SMm;
                g  = Gf' * eV_hat;

            else
                Of  = reshape(RVM(fbin,:,:), [V,M]);
                SMm = reshape(H_sm_model(fbin,:,:), [M,K]);
                SVm = reshape(H_sv_model(fbin,:,:), [V,K]);
                y_f = reshape(Yf(fbin,:), [K,1]);

                eV_hat = Of * (eMf - SMm*y_f) + SVm*y_f;
                g = SVm' * eV_hat;
            end

            xpow = sum(abs(Xvec).^2);
            eta = cfg.mu / (xpow + cfg.pwr_floor);

            for k = 1:K
                for j = 1:J
                    Wf(fbin,k,j) = (1-cfg.leak) * Wf(fbin,k,j) ...
                                 - eta * conj(Xvec(j)) * g(k);
                end
            end
        end
    end

    Dv_tf_on(:,1:tfrm_start-1,:) = Dv_tf_off(:,1:tfrm_start-1,:);
    Dv_tf_on(:,tfrm_end+1:end,:) = Dv_tf_off(:,tfrm_end+1:end,:);
    De_tf_on(:,1:tfrm_start-1,:) = De_tf_off(:,1:tfrm_start-1,:);
    De_tf_on(:,tfrm_end+1:end,:) = De_tf_off(:,tfrm_end+1:end,:);

    eval_off = zeros(Ns, P_eval);
    eval_on  = zeros(Ns, P_eval);

    for p = 1:P_eval
        eval_off(:,p) = force_len(istft(De_tf_off(:,:,p), fs, ...
            'Window', win, ...
            'OverlapLength', wlen-hop, ...
            'FFTLength', nfft, ...
            'FrequencyRange', 'onesided'), Ns);

        eval_on(:,p) = force_len(istft(De_tf_on(:,:,p), fs, ...
            'Window', win, ...
            'OverlapLength', wlen-hop, ...
            'FFTLength', nfft, ...
            'FrequencyRange', 'onesided'), Ns);
    end

    out.eval_off = eval_off;
    out.eval_on  = eval_on;
end