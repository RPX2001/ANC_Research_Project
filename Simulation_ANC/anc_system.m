% Helper Functions

function y = make_exact_length(x, N)
    if length(x) > N
        y = x(1:N);
    else
        % loop if shorter (better than zero-padding for noise)
        y = repmat(x, ceil(N/length(x)), 1);
        y = y(1:N);
    end
end
%% Parameters

fs_primary = 8000;
T_duration = 60;
Ns_target = fs_primary * T_duration;

Ns = size(src,1);                 % Size of the total samples
K  = size(src,2);                 % #reference channels (primary sources 2)
M  = size(h_pm,2);                % #monitoring mics (8)
Lspk = size(h_sm,3);              % #secondary speakers (2)
F = nfft/2 + 1;                   % frequency bins


% normalizing coefficients
ref_gain = 1 ./ (rms(src,1) + 1e-12);          % [1 x K]
mon_gain = 1 ./ (rms(mon_sig,1) + 1e-12);      % [1 x M]

% dimension change for ReTM
RVM = squeeze(retm_est(:,1,:,:));         % -> [F x V x M]
V = size(RVM,2);

%% Build frequency responses of primary and secondary paths
X = cell(K,1);
for k = 1:K
    X{k} = stft(src(:,k), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided');  % [F x nFrames]
end
nFrames = size(X{1},2);

% ---- FFT of paths
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

%% monitoring mic signal from primary to monitoring mics

% Primary contribution at monitoring mics: Dm_tf
Dm_tf = zeros(F, nFrames, M);
for k = 1:K
    for m = 1:M
        Dm_tf(:,:,m) = Dm_tf(:,:,m) + Hpm(:,m,k) .* X{k};
    end
end

% Primary contribution at evaluation mics (ANC OFF baseline): Deval_tf_off
Deval_tf_off = zeros(F, nFrames, Neval);
for k = 1:K
    for p = 1:Neval
        Deval_tf_off(:,:,p) = Deval_tf_off(:,:,p) + Hpee(:,p,k) .* X{k};
    end
end

%% ANC parameters

% Adaptive filter W(f,l,k)
W = zeros(F, Lspk, K);

% Store virtual error before/after + evaluation mics after
EV0_tf = zeros(F, nFrames, V);
EV_tf  = zeros(F, nFrames, V);
Deval_tf_on = zeros(F, nFrames, Neval);

% FxLMS stability settings (full-band needs smaller mu)
mu = 1e-3;              % start here for full-band
pwr_floor = 1e-3;       % prevents huge steps

adapt_bins = 2:(F-1);   % cover all the frequency range for the anc cancellation

%% Main ANC Loop

for tfrm = 1:nFrames
    for fbin = 1:F

        % Xvec [K x 1]
        Xvec = zeros(K,1);
        for k = 1:K
            Xvec(k) = X{k}(fbin,tfrm);
        end

        % Speaker outputs Y [L x 1]
        Y = zeros(Lspk,1);
        for l = 1:Lspk
            Y(l) = squeeze(W(fbin,l,:)).' * Xvec;
        end

        % Monitoring residual eM = dM + S_M Y
        dMf = squeeze(Dm_tf(fbin,tfrm,:));     % [M x 1]
        SMf = squeeze(Hsm(fbin,:,:));          % [M x L]
        eMf = dMf + SMf * Y;                   % [M x 1]

        % Monitoring normalization BEFORE ReTM (must match tuning)
        % Apply per-channel gain (vector) consistently
        eMfN = eMf .* mon_gain(:);

        % Virtual error (before) and (after)
        Rf = squeeze(RVM(fbin,:,:));           % [V x M]
        EV0_tf(fbin,tfrm,:) = Rf * (dMf .* mon_gain(:)); % virtual error without secondary speakers
        eVf = Rf * eMfN;
        EV_tf(fbin,tfrm,:) = eVf;              % estimated error signal

        % Evaluation mics (ANC ON): primary + secondary
        % e_eval = P_eval*X + S_eval*Y
        for p = 1:Neval
            primPart = 0;
            for k = 1:K
                primPart = primPart + Hpee(fbin,p,k) * Xvec(k);
            end
            secPart = squeeze(Hsee(fbin,p,:)).' * Y;
            Deval_tf_on(fbin,tfrm,p) = primPart + secPart;
        end

        % Weight update only for adapt bins
        if fbin >= adapt_bins(1) && fbin <= adapt_bins(end)

            % Gradient: g = S_M^H * R^H * eV
            g = (SMf') * (Rf') * eVf;          % [L x 1]

            % NLMS normalization
            pwr = real(Xvec' * Xvec);
            pwr = max(pwr, pwr_floor);

            for l = 1:Lspk
                for k = 1:K
                    W(fbin,l,k) = W(fbin,l,k) - (mu/pwr) * conj(Xvec(k)) * g(l);
                end
            end
        end

    end
end

%% ISTFT

ev_before = zeros(Ns, V);
ev_after  = zeros(Ns, V);

for v = 1:V
    ev_before(:,v) = istft(EV0_tf(:,:,v), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided');
    ev_after(:,v)  = istft(EV_tf(:,:,v), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided');
end

% ---- ISTFT: evaluation mics (ANC OFF / ON)
eval_off = zeros(Ns, Neval);
eval_on  = zeros(Ns, Neval);

for p = 1:Neval
    eval_off(:,p) = istft(Deval_tf_off(:,:,p), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided');
    eval_on(:,p)  = istft(Deval_tf_on(:,:,p), fs, ...
        'Window', win, 'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, 'FrequencyRange', 'onesided');
end

% Fix lengths
ev_before = ev_before(1:Ns,:);
ev_after  = ev_after(1:Ns,:);
eval_off  = eval_off(1:Ns,:);
eval_on   = eval_on(1:Ns,:);
%% plot

t = (0:Ns-1)'/fs;
pPlot = 1;

figure;
plot(t, eval_off(:,pPlot)); hold on;
plot(t, eval_on(:,pPlot));
grid on; xlabel('Time (s)'); ylabel('Amplitude');
legend('ANC OFF (eval mic)','ANC ON (eval mic)');
title(sprintf('Evaluation mic %d: time domain', pPlot));

% PSD plot and noise reduction
nwel=4096; nover=round(0.75*nwel); nfftW=4096;
trim=2*wlen;

xb = eval_off(trim:end-trim, pPlot);
xa = eval_on( trim:end-trim, pPlot);

[Pb,fpsd] = pwelch(xb, hann(nwel,'periodic'), nover, nfftW, fs);
[Pa,~]    = pwelch(xa, hann(nwel,'periodic'), nover, nfftW, fs);

Pb_dB = 10*log10(max(Pb,1e-20));
Pa_dB = 10*log10(max(Pa,1e-20));
NR = Pb_dB - Pa_dB;

% mask bins with no energy
thr = max(Pb_dB) - 40;
mask = Pb_dB > thr;
NR(~mask) = NaN;
NR = min(max(NR,-20),40);

figure;
plot(fpsd, Pb_dB); hold on; plot(fpsd, Pa_dB);
grid on; xlim([0 fs/2]);
xlabel('Hz'); ylabel('PSD (dB/Hz)');
legend('OFF','ON');
title(sprintf('Evaluation mic %d: PSD', pPlot));

figure;
plot(fpsd, NR, 'LineWidth', 1.2);
grid on; xlim([0 fs/2]);
xlabel('Hz'); ylabel('Noise reduction (dB)');
title(sprintf('Evaluation mic %d: Noise reduction vs frequency', pPlot));