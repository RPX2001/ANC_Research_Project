clearvars;
close all;
clc;

fprintf('=== Scenario 3: RM vs ReTM with secondary-path change ===\n');

function y = make_exact_length(x, N)
    if length(x) > N
        y = x(1:N);
    else
        % loop if shorter (better than zero-padding for noise)
        y = repmat(x, ceil(N/length(x)), 1);
        y = y(1:N);
    end
end

addpath(genpath('RIR_Generator'));
addpath(genpath('ReTM'));

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

num_pri_src = 2;
num_sec_src = 2;
num_mon_mics = 8;
num_err_mics = 4;
num_eval_err_mics = 2;

Ttune = 100;
Ns = fs*Ttune;
T = 100;
Ns_target = fs*T;

wlen = 1024;
hop  = 256;
nfft = 2048;
win  = hann(wlen,'periodic');
F    = nfft/2 + 1;

%% ============================================================
% Geometry
%% ============================================================
mon_mics = [0.15  0.15 -0.15;
            0.15 -0.15 -0.15;
           -0.15  0.15 -0.15;
           -0.15 -0.15 -0.15;
            0     0.212 0.15;
            0    -0.212 0.15;
            0.212 0     0.15;
           -0.212 0     0.15];

prim_sources = [0.6  0.8  1.00;
                2.17 1.15 0.05];

error_mics = [0.05  0.10 0;
              0.05 -0.10 0;
             -0.05  0.10 0;
             -0.05 -0.10 0];

eval_mic = [0  0.10 0;
            0 -0.10 0];

% Scenario 3: secondary-path change
sec_sources_tune = [0.00  0.35 0.00;
                    0.00 -0.35 0.00];

% Make the control-stage secondary path DIFFER noticeably
sec_sources_ctrl = [ 0  0.30 0.00;
                    0 -0.30 0.00];

%% ============================================================
% Source signals
%% ============================================================
fprintf('Loading source signals...\n');

[x1, fs1] = audioread('Sound_sources/buccaneer1.wav');
[x2, fs2] = audioread('Sound_sources/buccaneer2.wav');

if size(x1,2) > 1, x1 = mean(x1,2); end
if size(x2,2) > 1, x2 = mean(x2,2); end

if fs1 ~= fs, x1 = resample(x1, fs, fs1); end
if fs2 ~= fs, x2 = resample(x2, fs, fs2); end

x1 = make_exact_length(x1, Ns_target);
x2 = make_exact_length(x2, Ns_target);

src = [x1, x2];

%% ============================================================
% Tuning-stage impulse responses
%% ============================================================
fprintf('Generating tuning-stage impulse responses...\n');

h_pm = zeros(n, num_mon_mics, num_pri_src);
h_pe = zeros(n, num_err_mics, num_pri_src);
h_pee = zeros(n, num_eval_err_mics, num_pri_src);

for k = 1:num_pri_src
    h_pm(:,:,k) = rir_generator(c, fs, mon_mics, prim_sources(k,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
    h_pe(:,:,k) = rir_generator(c, fs, error_mics, prim_sources(k,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
    h_pee(:,:,k) = rir_generator(c, fs, eval_mic, prim_sources(k,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
end

h_sm_tune = zeros(n, num_mon_mics, num_sec_src);
h_se_tune = zeros(n, num_err_mics, num_sec_src);
h_see_tune = zeros(n, num_eval_err_mics, num_sec_src);

for l = 1:num_sec_src
    h_sm_tune(:,:,l) = rir_generator(c, fs, mon_mics, sec_sources_tune(l,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
    h_se_tune(:,:,l) = rir_generator(c, fs, error_mics, sec_sources_tune(l,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
    h_see_tune(:,:,l) = rir_generator(c, fs, eval_mic, sec_sources_tune(l,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
end

%% ============================================================
% Control-stage CHANGED secondary paths
%% ============================================================
fprintf('Generating control-stage changed secondary paths...\n');

h_sm_ctrl = zeros(n, num_mon_mics, num_sec_src);
h_se_ctrl = zeros(n, num_err_mics, num_sec_src);
h_see_ctrl = zeros(n, num_eval_err_mics, num_sec_src);

for l = 1:num_sec_src
    h_sm_ctrl(:,:,l) = rir_generator(c, fs, mon_mics, sec_sources_ctrl(l,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
    h_se_ctrl(:,:,l) = rir_generator(c, fs, error_mics, sec_sources_ctrl(l,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
    h_see_ctrl(:,:,l) = rir_generator(c, fs, eval_mic, sec_sources_ctrl(l,:), ...
        Lroom, beta, n, mtype, order, dim, orientation, hp_filter).';
end

%% ============================================================
% Tuning-stage monitoring/error signals
%% ============================================================
fprintf('Generating tuning-stage monitoring/error signals...\n');

sec_noise = randn(Ns, num_sec_src);
for l = 1:num_sec_src
    sec_noise(:,l) = bandpass(sec_noise(:,l), [20 600], fs);
end

X = cell(num_pri_src,1);
U = cell(num_sec_src,1);

for k = 1:num_pri_src
    X{k} = stft(src(:,k), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided');
end

for l = 1:num_sec_src
    U{l} = stft(sec_noise(:,l), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided');
end

nFrames = size(X{1},2);

H_pm = rir_to_H(h_pm, nfft);
H_pe = rir_to_H(h_pe, nfft);
H_sm_tune = rir_to_H(h_sm_tune, nfft);
H_se_tune = rir_to_H(h_se_tune, nfft);

EM_tf = zeros(F, nFrames, num_mon_mics);
EE_tf = zeros(F, nFrames, num_err_mics);

for k = 1:num_pri_src
    for m = 1:num_mon_mics
        EM_tf(:,:,m) = EM_tf(:,:,m) + H_pm(:,m,k) .* X{k};
    end
    for e = 1:num_err_mics
        EE_tf(:,:,e) = EE_tf(:,:,e) + H_pe(:,e,k) .* X{k};
    end
end

for l = 1:num_sec_src
    for m = 1:num_mon_mics
        EM_tf(:,:,m) = EM_tf(:,:,m) + H_sm_tune(:,m,l) .* U{l};
    end
    for e = 1:num_err_mics
        EE_tf(:,:,e) = EE_tf(:,:,e) + H_se_tune(:,e,l) .* U{l};
    end
end

mon_sig = zeros(Ns, num_mon_mics);
err_sig = zeros(Ns, num_err_mics);

for m = 1:num_mon_mics
    mon_sig(:,m) = force_len(istft(EM_tf(:,:,m), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided'), Ns);
end

for e = 1:num_err_mics
    err_sig(:,e) = force_len(istft(EE_tf(:,:,e), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided'), Ns);
end

%% ============================================================
% Estimate ReTM and RM
%% ============================================================
fprintf('Estimating ReTM...\n');
[retm_est, ~] = retm_estimate(err_sig, mon_sig, ...
    'fs', fs, 'wlen', wlen, 'nfft', nfft, 'hop', hop);

fprintf('Estimating RM...\n');
[rm_est, ~] = rm_estimate(err_sig, mon_sig, ...
    'fs', fs, 'wlen', wlen, 'nfft', nfft, 'hop', hop, 'reg', 1e-6);

%% ============================================================
% Scenario 3 control settings
%% ============================================================
cfg = struct();
cfg.t_start = 30;
cfg.t_end   = 70;
cfg.skipFrames = 60;
cfg.mu = 8e-4;
cfg.pwr_floor = 1e-6;
cfg.leak = 1e-5;
cfg.f_lo = 20;
cfg.f_hi = 600;
cfg.pPlot = 1;

%% ============================================================
% Run ReTM
%% ============================================================
fprintf('\nRunning Scenario 3 with ReTM...\n');
outReTM = run_scenario3_control( ...
    'ReTM', retm_est, ...
    h_pm, h_pee, ...
    h_sm_ctrl, h_see_ctrl, ...
    H_sm_tune, H_se_tune, ...
    src, fs, wlen, hop, nfft, win, cfg);

%% ============================================================
% Run RM
%% ============================================================
fprintf('\nRunning Scenario 3 with RM...\n');
outRM = run_scenario3_control( ...
    'RM', rm_est, ...
    h_pm, h_pee, ...
    h_sm_ctrl, h_see_ctrl, ...
    H_sm_tune, H_se_tune, ...
    src, fs, wlen, hop, nfft, win, cfg);

%% ============================================================
% PSD / NR plots
%% ============================================================
p = cfg.pPlot;
trim = 2*wlen;
nwel = 4096;
nover = round(0.75*nwel);
nfftW = 4096;

xb = outReTM.eval_off(trim:end-trim,p);
xa_retm = outReTM.eval_on(trim:end-trim,p);
xa_rm   = outRM.eval_on(trim:end-trim,p);

[Pb, f] = pwelch(xb, hann(nwel,'periodic'), nover, nfftW, fs);
[Pa_retm, ~] = pwelch(xa_retm, hann(nwel,'periodic'), nover, nfftW, fs);
[Pa_rm, ~]   = pwelch(xa_rm, hann(nwel,'periodic'), nover, nfftW, fs);

Pb_dB      = 10*log10(max(Pb, 1e-20));
Pa_retm_dB = 10*log10(max(Pa_retm, 1e-20));
Pa_rm_dB   = 10*log10(max(Pa_rm, 1e-20));

NR_retm = Pb_dB - Pa_retm_dB;
NR_rm   = Pb_dB - Pa_rm_dB;

figure('Name','Scenario 3 PSD');
plot(f, Pa_retm_dB, 'k', 'LineWidth', 1.5); hold on;
plot(f, Pa_rm_dB, '--m', 'LineWidth', 1.5);
plot(f, Pb_dB, ':b', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Power / PSD (dB/Hz)');
legend('ReTM','RM','No ANC','Location','best');
title('Scenario 3: With secondary-path change');
xlim([0 600]);

figure('Name','Scenario 3 Noise Reduction');
plot(f, NR_retm, 'k', 'LineWidth', 1.5); hold on;
plot(f, NR_rm, '--m', 'LineWidth', 1.5);
yline(0,'k:');
grid on;
xlabel('Frequency (Hz)');
ylabel('Noise Reduction (dB)');
legend('ReTM','RM','Location','best');
title('Scenario 3: Noise reduction');
xlim([0 600]);

fprintf('\nAvg NR 50-400 Hz: ReTM = %.2f dB, RM = %.2f dB\n', ...
    nanmean(NR_retm(f>=50 & f<=400)), nanmean(NR_rm(f>=50 & f<=400)));

%% ============================================================
% Local helper functions
%% ============================================================

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

function y = force_len(x, N)
x = x(:);
if numel(x) >= N
    y = x(1:N);
else
    y = [x; zeros(N-numel(x),1)];
end
end

function out = run_scenario3_control(methodName, map_est, ...
    h_pm, h_pee, h_sm_ctrl, h_see_ctrl, ...
    H_sm_tune, H_se_tune, ...
    src_full, fs, wlen, hop, nfft, win, cfg)

F = nfft/2 + 1;
t_start = cfg.t_start;
t_end   = cfg.t_end;
n_start = floor(t_start*fs) + 1;
n_end   = floor(t_end*fs);

src = src_full(n_start:n_end,:);
Ns = size(src,1);
K = size(src,2);
M = size(h_pm,2);
Lspk = size(h_sm_ctrl,3);
Neval = size(h_pee,2);
V = size(map_est,3);

RVM = reshape(map_est(:,1,:,:), [F, V, M]);

% STFT of references
X = cell(K,1);
for k = 1:K
    X{k} = stft(src(:,k), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided');
end
nFrames = size(X{1},2);

% Actual control-stage paths
Hpm = rir_to_H(h_pm, nfft);
Hsm_ctrl = rir_to_H(h_sm_ctrl, nfft);
Hpee = rir_to_H(h_pee, nfft);
Hsee_ctrl = rir_to_H(h_see_ctrl, nfft);

% Primary fields
Dm_tf = zeros(F, nFrames, M);
Deval_tf_off = zeros(F, nFrames, Neval);

for k = 1:K
    for m = 1:M
        Dm_tf(:,:,m) = Dm_tf(:,:,m) + Hpm(:,m,k) .* X{k};
    end
    for p = 1:Neval
        Deval_tf_off(:,:,p) = Deval_tf_off(:,:,p) + Hpee(:,p,k) .* X{k};
    end
end

% Controller state
Wf = zeros(F, Lspk, K);
EV_tf = zeros(F, nFrames, V);
Deval_tf_on = zeros(F, nFrames, Neval);

Xmat = zeros(F, nFrames, K);
for k = 1:K
    Xmat(:,:,k) = X{k};
end

f_axis = (0:F-1)' * fs/nfft;
adapt_bins = find(f_axis >= cfg.f_lo & f_axis <= cfg.f_hi);

skipFrames = cfg.skipFrames;
tfrm_start = 1 + skipFrames;
tfrm_end   = nFrames - skipFrames;

for tfrm = tfrm_start:tfrm_end

    % Speaker outputs
    Yf = zeros(F, Lspk);
    for l = 1:Lspk
        acc = zeros(F,1);
        for k = 1:K
            acc = acc + Wf(:,l,k) .* Xmat(:,tfrm,k);
        end
        Yf(:,l) = acc;
    end

    % Actual plant uses CHANGED secondary path
    eM_tf_frame = squeeze(Dm_tf(:,tfrm,:));
    for l = 1:Lspk
        eM_tf_frame = eM_tf_frame + squeeze(Hsm_ctrl(:,:,l)) .* Yf(:,l);
    end

    % Virtual error estimate
    for fbin = 1:F
        eMf = reshape(eM_tf_frame(fbin,:), [M,1]);

        if strcmpi(methodName, 'ReTM')
            Rf = reshape(RVM(fbin,:,:), [V,M]);
            eVf = Rf * eMf;
        else
            Of = reshape(RVM(fbin,:,:), [V,M]);
            SM_t = reshape(H_sm_tune(fbin,:,:), [M,Lspk]);
            SV_t = reshape(H_se_tune(fbin,:,:), [V,Lspk]);
            y_f = reshape(Yf(fbin,:), [Lspk,1]);

            eVf = Of * (eMf - SM_t*y_f) + SV_t*y_f;
        end

        EV_tf(fbin,tfrm,:) = eVf;
    end

    % Evaluation microphones use CHANGED path
    Deval_frame = squeeze(Deval_tf_off(:,tfrm,:));
    for l = 1:Lspk
        Deval_frame = Deval_frame + squeeze(Hsee_ctrl(:,:,l)) .* Yf(:,l);
    end
    Deval_tf_on(:,tfrm,:) = Deval_frame;

    % Controller update uses OLD / TUNING path model
    for ib = 1:numel(adapt_bins)
        fbin = adapt_bins(ib);

        Xvec = reshape(Xmat(fbin,tfrm,:), [K,1]);
        eVf = reshape(EV_tf(fbin,tfrm,:), [V,1]);

        SM_model = reshape(H_sm_tune(fbin,:,:), [M,Lspk]);   % model used in controller

        tmp = eVf_to_monitor(methodName, fbin, eVf, RVM, H_sm_tune, H_se_tune, Yf, M, V, Lspk);
        g = (SM_model') * tmp;

        xpow = sum(abs(Xvec).^2);
        eta = cfg.mu / (xpow + cfg.pwr_floor);

        for l = 1:Lspk
            for k = 1:K
                Wf(fbin,l,k) = (1-cfg.leak) * Wf(fbin,l,k) - eta * conj(Xvec(k)) * g(l);
            end
        end
    end
end

Deval_tf_on(:,1:tfrm_start-1,:) = Deval_tf_off(:,1:tfrm_start-1,:);
Deval_tf_on(:,tfrm_end+1:end,:) = Deval_tf_off(:,tfrm_end+1:end,:);

eval_off = zeros(Ns,Neval);
eval_on  = zeros(Ns,Neval);
for p = 1:Neval
    eval_off(:,p) = force_len(istft(Deval_tf_off(:,:,p), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided'), Ns);
    eval_on(:,p)  = force_len(istft(Deval_tf_on(:,:,p), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided'), Ns);
end

out.eval_off = eval_off;
out.eval_on  = eval_on;
out.t = ((n_start-1)+(0:Ns-1))'/fs;
end

function tmp = eVf_to_monitor(methodName, fbin, eVf, RVM, H_sm_tune, H_se_tune, Yf, M, V, Lspk)
    if strcmpi(methodName, 'ReTM')
        Rf = reshape(RVM(fbin,:,:), [V, M]);   % [V x M]
        tmp = Rf' * eVf;                       % [M x 1]
    else
        Of = reshape(RVM(fbin,:,:), [V, M]);   % [V x M]
        tmp = Of' * eVf;                       % [M x 1]
    end
    
    tmp = reshape(tmp, [M, 1]);
end