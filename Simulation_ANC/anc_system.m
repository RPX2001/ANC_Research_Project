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

% ================= CONTROL STAGE TIME WINDOW =================
t_start = 20;              % seconds
t_end   = 40;              % seconds
n_start = floor(t_start*fs) + 1;
n_end   = floor(t_end*fs);

% Keep a copy of full source if you need it later
src_full = src;

% Crop to 10–50 s for control stage ONLY
src = src_full(n_start:n_end, :);

Ns = size(src,1);                 % Size of the total samples
K  = size(src,2);                 % #reference channels (primary sources 2)
M  = size(h_pm,2);                % #monitoring mics (8)
Lspk = size(h_sm,3);              % #secondary speakers (2)
F = nfft/2 + 1;                   % frequency bins
Neval = size(h_pee,2);            % # of evaluation mics - 2
W_nums = 64;                      % filter taps

% history buffer for reference signals (time domain)
x_hist = zeros(W_nums, K);          % last 64 samples per reference channel

% normalizing coefficients
ref_gain = ones(1,K);      % no reference normalization
mon_gain = ones(1,M);      % no monitoring normalization

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
W_nums = 64;                      % Filter coefficients
w_td = zeros(Lspk, K, W_nums);    % time-domain FIR for each (speaker,ref)
Wf   = zeros(F, Lspk, K);         % frequency response of each FIR tap set

hopS = hop;                       % samples per frame hop (256)


% Store virtual error before/after + evaluation mics after
EV0_tf = zeros(F, nFrames, V);
EV_tf  = zeros(F, nFrames, V);
Deval_tf_on = zeros(F, nFrames, Neval);

% FxLMS stability settings (full-band needs smaller mu)
mu = 5e-4;              % start here for full-band
pwr_floor = 1e-3;       % prevents huge steps

adapt_bins = 2:(F-1);   % cover all the frequency range for the anc cancellation

%% Main ANC Loop

% Stack references into one array for easy access:
% Xmat: [F x nFrames x K]
Xmat = zeros(F, nFrames, K);
for k = 1:K
    Xmat(:,:,k) = X{k};
end


% Skip first/last frames (optional safety)
skipFrames = 60; % try 5..20
tfrm_start = 1 + skipFrames;
tfrm_end   = nFrames - skipFrames;

for tfrm = tfrm_start:tfrm_end

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

    % --- STFT-domain FxLMS update (paper-style) ---
    for fbin = adapt_bins

        % reference vector x(f,t) = [X1; X2; ... XK]
        Xvec = squeeze(Xmat(fbin,tfrm,:));        % [K x 1]

        % compute normalization power
        pwr = sum(abs(Xvec).^2);
        pwr = max(pwr, pwr_floor);

        % get error vector at virtual mics
        eVf = squeeze(EV_tf(fbin,tfrm,:));        % [V x 1]

        % Secondary path to monitoring mics (M x L)
        SMf = squeeze(Hsm(fbin,:,:));             % [M x L]

        % ReTM mapping (V x M)
        Rf = squeeze(RVM(fbin,:,:));              % [V x M]

        % gradient term per speaker:
        % g = S_M^H * R^H * eV
        g = (SMf') * (Rf') * eVf;                 % [L x 1]

        % Update weights
        for l = 1:Lspk
            for k = 1:K
                Wf(fbin,l,k) = Wf(fbin,l,k) - (mu/pwr) * conj(Xvec(k)) * g(l);
            end
        end
    end

end

% Fill skipped frames to avoid ISTFT artifacts
EV_tf(:,1:tfrm_start-1,:) = EV0_tf(:,1:tfrm_start-1,:);
EV_tf(:,tfrm_end+1:end,:) = EV0_tf(:,tfrm_end+1:end,:);
Deval_tf_on(:,1:tfrm_start-1,:) = Deval_tf_off(:,1:tfrm_start-1,:);
Deval_tf_on(:,tfrm_end+1:end,:) = Deval_tf_off(:,tfrm_end+1:end,:);




%% ISTFT

force_len = @(x,N) (length(x)>=N)*x(1:N) + (length(x)<N)*[x(:); zeros(N-length(x),1)];

% ISTFT virtual error
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

% ISTFT evaluation mics
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

%% plot

t = ((n_start-1) + (0:Ns-1))' / fs;   % now t runs from 10 to 50 seconds
pPlot = 1;

figure;
plot(t, eval_off(:,pPlot)); hold on;
plot(t, eval_on(:,pPlot));
grid on; xlabel('Time (s)'); ylabel('Amplitude');
legend('ANC OFF (eval mic)','ANC ON (eval mic)');
title(sprintf('Evaluation mic %d: time domain', pPlot));
xlim([25 35]);

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
thr = max(Pb_dB) - 100; % 40
mask = Pb_dB > thr;
NR(~mask) = NaN;
NR = min(max(NR,-30),60); % 20

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