% Example of ReTF Denoise
% Remove a unmoving continuous noise source from a 2-channel recording.
%
% Inputs:
%   target_recording.wav  (or the first channel of 2ch recording).
%   reference_recording.wav  (or the secondd channel of 2ch recording).
% 
% Outputs:
%   retf - Estimated ReTF.
%   denoised_sig - ReTF Denoised target signal.
%
% Other m-files required: +shaasp functions.
% Subfunctions: none
% MAT-files required: none
%
% See also: example_retm_denoise
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 02-Oct-2024
% Last revision: 09-Oct-2024

% --- Step 1) User Inputs ---------------

% Load target signal / channel.
[tgt_sig, tgt_fs] = audioread('test.wav');
tgt_sig = tgt_sig(:,1);

% Load reference signal / channel.
[ref_sig, ref_fs] = audioread('test.wav');
ref_sig = ref_sig(:,2);

% ---------------------------------------

% Pre-process settings.
fs = 48000;

% STFT settings.
wlen = fs / 1;  % Window length should be longer than RT60.
nfft = 2^nextpow2(2 * wlen);  % NFFT size should be greater than 2 * wlen.
hop = wlen / 4;  % Step / hop size should be <= wlen.
ola_method = 'ola';  % ReTF is filtering so need STFT OLA processing.
wind = shaasp.cola_window(wlen, hop, ola_method);  % COLA hann window.

% Pre-process the signals, mono, resample, normalize, match sizes.
tgt_sig = tgt_sig(:,1);
ref_sig = ref_sig(:,1);

tgt_sig = resample(tgt_sig, fs, tgt_fs);
ref_sig = resample(ref_sig, fs, ref_fs);

tgt_sig = tgt_sig ./ max(abs(tgt_sig(:)));
ref_sig = ref_sig ./ max(abs(ref_sig(:)));

tgt_nsamples = length(tgt_sig);
ref_nsamples = length(ref_sig);

if (ref_nsamples < tgt_nsamples)
    ref_sig = [ref_sig; zeros(tgt_nsamples - ref_nsamples, 1)];
elseif (ref_nsamples > tgt_nsamples)
    ref_sig = ref_sig(1:tgt_nsamples, :);
end

% Step 2) Short Time Fourier Transform (STFT) spectrum of the signals.
tgt_spec = shaasp.cola_stft(tgt_sig, nfft, hop, wind);  % [f,t,1]
ref_spec = shaasp.cola_stft(ref_sig, nfft, hop, wind);  % [f,t,1]

% Step 3) Cross Power Spectral Density (PSD). 
psd_tgt_ref = shaasp.spec_to_crossspec(tgt_spec, ref_spec);  % [f,t,1,1]
psd_ref_ref = shaasp.spec_to_crossspec(ref_spec);  % auto-correlation.

% Step 4) Spatial correlation (time averaged PSDs).
averaging_method = 'all';
% % averaging_method = 'shorttime';

switch averaging_method
    case 'all'
        % Average over all time frames for time-invariant ReTF [f,1,1].
        corr_tgt_ref = shaasp.spec_time_averaging(psd_tgt_ref, [], "all");
        corr_ref_ref = shaasp.spec_time_averaging(psd_ref_ref, [], "all");
    
    case 'shorttime'
        % Average over time window for shorttime ReTF [f,t,1]. 
        ave_time_seconds = 5;
        ave_time_nframes = floor(ave_time_seconds * fs / hop);
        corr_tgt_ref = shaasp.spec_time_averaging(psd_tgt_ref, ave_time_nframes, "closest");
        corr_ref_ref = shaasp.spec_time_averaging(psd_ref_ref, ave_time_nframes, "closest");

    otherwise
        error('Typo somewhere :)');
end

% Step 5) Relative Transfer Function (ReTF)
retf = corr_tgt_ref ./ corr_ref_ref;

% Step 6) Estimate noise component spectrum from the reference.
est_spec = retf .* ref_spec;

% ISTFT to get noise component signal in target signal.
est_sig = shaasp.cola_istft(est_spec, nfft, hop);
est_sig(isnan(est_sig)) = 0;
est_sig = est_sig(1 : tgt_nsamples);

% Step 7) Denoise target signal by subtracting the estimated noise component.
denoised_sig = tgt_sig - est_sig;

return