clear all;
close all;
%% Helper functions

function y = make_exact_length(x, N)
    if length(x) > N
        y = x(1:N);
    else
        % loop if shorter (better than zero-padding for noise)
        y = repmat(x, ceil(N/length(x)), 1);
        y = y(1:N);
    end
end


%% RIR Generation

c = 340;                    % Sound velocity (m/s)
fs = 8000;                  % Sample frequency (samples/s)
L = [8 7 4];                % Room dimensions [x y z] (m)
beta = 0.4;                 % Reverberation time (s) - adjust as needed
n = 4096;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                 % -1 equals maximum reflection order
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

num_pri_src = 2;            % number of primary sources
num_sec_src = 2;            % number of secondary speakers
num_mon_mics = 8;           % number of monitoring mics
num_err_mics = 4;           % number of error mics
num_eval_err_mics = 2;      % number of evaluation error mics

%fs = 8000;                  % secondary source target sampling rate
Ttune = 60;                 % Duration of the signal      
Ns = fs*Ttune;              % total samples in secondary source                                        

%fs_target = 8000;           % primary source target sampling rate
T = 60;                     % primary source signal duration seconds
Ns_target = fs * T;         % primary source total samples

wlen = 1024;
hop  = 256;
nfft = 2048;
win  = hann(wlen,'periodic');
F    = nfft/2 + 1;

%% Monitoring microphones positions

mon_mics = [0.15  0.15 -0.15;           % Microphone 1
            0.15 -0.15 -0.15;           % Microphone 2
           -0.15  0.15 -0.15;           % Microphone 3
           -0.15 -0.15 -0.15;           % Microphone 4
            0     0.212 0.15;           % mic 5
            0    -0.212 0.15;           % mic 6
            0.212 0     0.15;           % mic 7
           -0.212 0     0.15];          % mic 8
%% Primary sources and secondary sources positions

prim_sources = [0.6  0.8  1;
                2.17 1.15 0.05];  % primary sources

sec_sources = [0  0.3 0;
               0 -0.3 0];        % Secondary sources

%% Error microphone positions

error_mics = [0.05  0.10    0;
              0.05 -0.10    0;
             -0.05  0.10    0;
             -0.05 -0.10    0];   % error microphones

eval_mic = [0   0.10   0;
            0  -0.10   0];        % ANC evaluation mics

%% Generate Impulese Responses

h_pm = zeros(n, num_mon_mics, num_pri_src);

for k = 1:num_pri_src
    h_pm(:,:,k) = rir_generator( ...
        c, fs, mon_mics, prim_sources(k,:), ...
        L, beta, n, mtype, order, dim, orientation, hp_filter).';
end

h_pe = zeros(n, num_err_mics, num_pri_src);

for k = 1:num_pri_src
    h_pe(:,:,k) = rir_generator( ...
        c, fs, error_mics, prim_sources(k,:), ...
        L, beta, n, mtype, order, dim, orientation, hp_filter).';
end

h_sm = zeros(n, num_mon_mics, num_sec_src);

for l = 1:num_sec_src
    h_sm(:,:,l) = rir_generator( ...
        c, fs, mon_mics, sec_sources(l,:), ...
        L, beta, n, mtype, order, dim, orientation, hp_filter).';
end

h_se = zeros(n, num_err_mics, num_sec_src);

for l = 1:num_sec_src
    h_se(:,:,l) = rir_generator( ...
        c, fs, error_mics, sec_sources(l,:), ...
        L, beta, n, mtype, order, dim, orientation, hp_filter).';
end

% Impulse response for evaluation error microphone from primary sources and
% secondary sources 

h_pee = zeros(n, num_eval_err_mics, num_pri_src);
h_see = zeros(n, num_eval_err_mics, num_sec_src);

for l = 1:num_sec_src
    h_see(:,:,l) = rir_generator( ...
        c, fs, eval_mic, sec_sources(l,:), ...
        L, beta, n, mtype, order, dim, orientation, hp_filter).';
end

for k = 1:num_pri_src
    h_pee(:,:,k) = rir_generator( ...
        c, fs, eval_mic, prim_sources(k,:), ...
        L, beta, n, mtype, order, dim, orientation, hp_filter).';
end


%% Generate Primary source signal -  Tuning Stage

% Load audio files
[x1, fs1] = audioread('factory1.wav');
[x2, fs2] = audioread('factory2.wav');

% Convert to mono if needed
if size(x1,2) > 1
    x1 = mean(x1,2);
end
if size(x2,2) > 1
    x2 = mean(x2,2);
end

% Make exactly 60 seconds
x1 = make_exact_length(x1, Ns_target);
x2 = make_exact_length(x2, Ns_target);

% Remove DC (important for ANC stability)
x1 = x1 - mean(x1);
x2 = x2 - mean(x2);

% Normalize power (do NOT normalize peak)
x1 = x1 / rms(x1);
x2 = x2 / rms(x2);

% Pack as primary sources
src = [x1, x2];    % size: [Ns_target x 2]

fprintf('Primary sources ready: %d samples @ %d Hz\n', size(src,1), fs);
%% Genrate secondary source signal - Tuning stage

sec_noise = randn(Ns, num_sec_src);

% Band-limit to control band
band = [1500 2400];
for l = 1:num_sec_src
    sec_noise(:,l) = bandpass(sec_noise(:,l), band, fs);
end

% Normalize (important)
sec_noise = sec_noise ./ (rms(sec_noise,1) + 1e-8);

% Scale: secondary noise should be weaker than primary
sec_gain = 0.3;                   
sec_noise = sec_gain * sec_noise;

%% Generate Monitoring microphone signals - Tuning Stage

% e_m = h_pm * x_t_f + h_sm * y_t_f

X = cell(num_pri_src,1);           % primary STFTs
U = cell(num_sec_src,1);           % secondary STFTs

for k = 1:num_pri_src
    X{k} = stft(src(:,k), fs, ...
        'Window', win, ...
        'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, ...
        'FrequencyRange','onesided');         % [F x nFrames]
end

for l = 1:num_sec_src
    U{l} = stft(sec_noise(:,l), fs, ...
        'Window', win, ...
        'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, ...
        'FrequencyRange','onesided');
end

nFrames = size(X{1},2);

% Convert each RIR to [Ch x n], then FFT to [F x Ch]

% Primary -> Monitoring: H_pm [F x M x K]
H_pm = zeros(F, num_mon_mics, num_pri_src);
for k = 1:num_pri_src
    Ht = squeeze(h_pm(:,:,k)).';    % [M × n]
    for m = 1:num_mon_mics
        H = fft(Ht(m,:), nfft);
        H_pm(:,m,k) = H(1:F).';
    end
end

% Primary -> Error: H_pe [F x E x K]
H_pe = zeros(F, num_err_mics, num_pri_src);
for k = 1:num_pri_src
    Ht = squeeze(h_pe(:,:,l)).';    % [E × n]
    for e = 1:num_err_mics
        H = fft(Ht(e,:), nfft);
        H_pe(:,e,k) = H(1:F).';
    end
end

% Secondary -> Monitoring: H_sm [F x M x L]
H_sm = zeros(F, num_mon_mics, num_sec_src);
for l = 1:num_sec_src
    Ht = squeeze(h_sm(:,:,k)).';    % [M × n]
    for m = 1:num_mon_mics
        H = fft(Ht(m,:), nfft);
        H_sm(:,m,l) = H(1:F).';
    end
end

% Secondary -> Error: H_se [F x E x L]
H_se = zeros(F, num_err_mics, num_sec_src);
for l = 1:num_sec_src
    Ht = squeeze(h_se(:,:,k)).';    % [M × n]
    for e = 1:num_err_mics
        H = fft(Ht(e,:), nfft);
        H_se(:,e,l) = H(1:F).';
    end
end

% Monitoring STFT: EM_tf [F x nFrames x M]
% Error STFT:      EE_tf [F x nFrames x E]

EM_tf = zeros(F, nFrames, num_mon_mics);
EE_tf = zeros(F, nFrames, num_err_mics);

% ---- primary contributions ----
for k = 1:num_pri_src
    Xk = X{k};  % [F x nFrames]
    for m = 1:num_mon_mics
        EM_tf(:,:,m) = EM_tf(:,:,m) + H_pm(:,m,k) .* Xk;
    end
    for e = 1:num_err_mics
        EE_tf(:,:,e) = EE_tf(:,:,e) + H_pe(:,e,k) .* Xk;
    end
end

% ---- secondary contributions ----
for l = 1:num_sec_src
    Ul = U{l};  % [F x nFrames]
    for m = 1:num_mon_mics
        EM_tf(:,:,m) = EM_tf(:,:,m) + H_sm(:,m,l) .* Ul;
    end
    for e = 1:num_err_mics
        EE_tf(:,:,e) = EE_tf(:,:,e) + H_se(:,e,l) .* Ul;
    end
end


% get the final monitoring and error signal

mon_sig = zeros(Ns, num_mon_mics);
err_sig = zeros(Ns, num_err_mics);

for m = 1:num_mon_mics
    mon_sig(:,m) = istft(EM_tf(:,:,m), fs, ...
        'Window', win, ...
        'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, ...
        'FrequencyRange','onesided');
end

for e = 1:num_err_mics
    err_sig(:,e) = istft(EE_tf(:,:,e), fs, ...
        'Window', win, ...
        'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, ...
        'FrequencyRange','onesided');
end

% Trim/force length 
mon_sig = mon_sig(1:Ns,:);
err_sig = err_sig(1:Ns,:);

disp(size(mon_sig));   
disp(size(err_sig));   

%% Add these two signal to ReTM estimator 

[retm_est, reim_est] = retm_estimate( ...
    err_sig, ...
    mon_sig, ...
    'fs', fs, ...
    'wlen', wlen, ...
    'nfft', nfft, ...
    'hop', hop);
%% freqency map


fbin = 100;          % frequency index
tfrm = 1;            % time/frame index (or 1 if it’s constant)

Re = squeeze(reim_est(fbin, tfrm, :, :));   % -> [nErr x nMon]


%% Recreate the error signal based on monitoring signal for varification

% STFT of monitoring signals
win = hann(wlen,'periodic');

R = cell(num_mon_mics,1);
for m = 1:num_mon_mics
    R{m} = stft(mon_sig(:,m), fs, ...
        'Window', win, ...
        'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, ...
        'FrequencyRange','onesided');   % [F x nFrames]
end
F = size(R{1},1);
nFrames = size(R{1},2);

% Apply ReTM in TF domain to estimate error STFT

Ehat_tf = zeros(F, nFrames, num_err_mics);

for e = 1:num_err_mics
    for m = 1:num_mon_mics
        H_em = squeeze(retm_est(:,1,e,m));   % [F x 1] (expected)
        H_em = H_em(:);
        Ehat_tf(:,:,e) = Ehat_tf(:,:,e) + H_em .* R{m};
    end
end

% ISTFT back

err_hat = zeros(size(err_sig));

for e = 1:num_err_mics
    err_hat(:,e) = istft(Ehat_tf(:,:,e), fs, ...
        'Window', win, ...
        'OverlapLength', wlen-hop, ...
        'FFTLength', nfft, ...
        'FrequencyRange','onesided');
end

% force same length (ISTFT can be off by a few samples)
Ns = size(err_sig,1);
err_hat = err_hat(1:Ns,:);

% plot error
ch = 2;
t = (0:Ns-1)/fs;

figure;
plot(t, err_sig(:,ch) - err_hat(:,ch));
grid on;
xlabel('Time (s)'); ylabel('Error');
title(sprintf('Estimation Error (True - Estimated), Mic %d', ch));

xlim([100 1000]);  


