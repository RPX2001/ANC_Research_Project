% Example of ReTM Denoise
% Remove multiple unmoving continuous noise sources from a multi-channel 
% recording.
%
% Inputs:
%   tgt_sig - target mics/channels signals to be denoise.
%   ref_sig - refernce mics/channels signals to be ReTM filtered.
%   tgt_training_sig - target mics/channels signals to learn ReTM from.
%   ref_training_sig - reference mics/channels signals to learn ReTM from. 
% 
% Outputs:
%   retf - Estimated ReTM.
%   denoised_sig - ReTM Denoised target signal.
%
% Other m-files required: +shaasp functions.
% Subfunctions: lastpagepinv(), lastpagemtimes()
% MAT-files required: none
%
% See also: example_retf_denoise
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 04-Oct-2024
% Last revision: 09-Oct-2024

% --- User Load Signal to Denoise -----------------------------------------
% Load recording of target and reference microphones to be denoised.
% This section creates the target and reference signals.
%
% Target signals (tgt_sig) are the microphone signals of the target group.
% The target signals will be denoised.
%
% References signals (ref_sig) are the microphone signals of the reference
% group. The reference signals will be filtered by the ReTM to estimate the
% noise at the target group.

% Load target signals / channels.
[tgt_sig, tgt_fs] = audioread('test2.wav');
tgt_sig = tgt_sig(:, (1:2:end));  % Allocate channels for target.

% Load reference signals / channels.
[ref_sig, ref_fs] = audioread('test2.wav');
ref_sig = ref_sig(:, (2:2:end));  % Allocate channels for reference.

% --- User Load Signal to Train ReTM --------------------------------------
% Denote / load recordings of target and reference for estimating ReTM.
% This recording should only contain active noise sources, and not contain
% any speech you want to enhance.

% Option 1) Learn the ReTM from the same recording that will be denoised.
% Assuming that the majority of the recording is just noise only.
retm_method = 'embedded';
tgt_training_sig = tgt_sig;
ref_training_sig = ref_sig;

% Option 2) Learn the ReTM from a different recording from the same setup.
% % retm_method = 'separate';   
% % [tgt_training_sig, tgt_fs] = audioread('.wav');
% % tgt_training_sig = tgt_training_sig(:, (1:2:end));  % Allocate training target.
% % [ref_training_sig, ref_fs] = audioread('.wav');
% % ref_training_sig = ref_training_sig(:, (2:2:end));  % Allocate training reference.

% -------------------------------------------------------------------------

% Pre-processing settings.
fs = 16000;  % Resample all signals to this.

% STFT settings.
wlen = fs / 2;  % Window length should be longer than RT60.
nfft = 2^nextpow2(2 * wlen);  % NFFT size should be greater than 2 * wlen.
hop = wlen / 4;  % Step / hop size should be <= wlen/2.
ola_method = 'ola';  % ReTF is filtering so need STFT OLA processing.
wind = shaasp.cola_window(wlen, hop, ola_method);  % COLA hann window.
single_sided = true;  % Single side to reduce computation.

% ReTM settings.
correlation_time_sec = 20;  % Seconds to average correlation over.
retm_smooth_time_sec = 20;  % Seconds to smooth ReTM over.

%% Pre-process signals.

% Resample to fs.
tgt_sig = resample(tgt_sig, fs, tgt_fs);
ref_sig = resample(ref_sig, fs, ref_fs);
tgt_training_sig = resample(tgt_training_sig, fs, tgt_fs);
ref_training_sig = resample(ref_training_sig, fs, tgt_fs);

% Normalize +-1 with equal gain on all recordings.
norm_term = max(abs([tgt_sig(:); ref_sig(:); tgt_training_sig(:); ref_training_sig(:)]));
tgt_sig = tgt_sig ./ norm_term;
ref_sig = ref_sig ./ norm_term;
tgt_training_sig = tgt_training_sig ./ norm_term;
ref_training_sig = ref_training_sig ./ norm_term;

% Trim / pad reference to same length as target.
ref_sig = trim_ref_to_tgt(ref_sig, tgt_sig);
ref_training_sig = trim_ref_to_tgt(ref_training_sig, tgt_training_sig);

%% Estimate ReTM.

% Short Time Fourier Transform (STFT) spectrum of training signals.
tgt_training_spec = shaasp.cola_stft(tgt_training_sig, nfft, hop, wind, single_sided);  % [f,t,a]
ref_training_spec = shaasp.cola_stft(ref_training_sig, nfft, hop, wind, single_sided);

% Cross-Power Spectral Density, 'a' = target, 'b' = reference.
cpsd_aa = shaasp.spec_to_crossspec(tgt_training_spec);  % auto-correlation. [f,t,a,a]
cpsd_ba = shaasp.spec_to_crossspec(ref_training_spec, tgt_training_spec);  % [f,t,b,a]

% Spatial correlation (time averaged CPSDs).
switch retm_method
    
    case 'embedded'
        % Estimate ReTM from averaging the mixture recording.
        ave_time_nframes = floor(correlation_time_sec * fs / hop);
        corr_aa = shaasp.spec_time_averaging(cpsd_aa, ave_time_nframes, "closest");
        corr_ba = shaasp.spec_time_averaging(cpsd_ba, ave_time_nframes, "closest");

    otherwise
        % Esimate ReTM as mean of all training CPSD frames.
        corr_aa = shaasp.spec_time_averaging(cpsd_aa, [], "all");
        corr_ba = shaasp.spec_time_averaging(cpsd_ba, [], "all");
end

% Estimate inverse of corr_ba. 
corr_ba_inv = lastpagepinv(corr_ba);  % [f,t,a,b]

% Estimate ReTM by ReTM = corr_aa * (corr_ba)^-1
retm_ab = lastpagemtimes(corr_aa, corr_ba_inv);  % [f,t,a,b]

%% Need to smooth the ReTM after matrix devision (not sure why).
switch retm_method
    case 'embedded'
        if ~isempty(retm_smooth_time_sec)
            retm_smooth_nframes = floor(retm_smooth_time_sec * fs / hop);
            retm_ab = shaasp.spec_time_averaging(retm_ab, retm_smooth_nframes, "closest");
        end

    otherwise
end

%% Denoise target signal using ReTM.
% Get stft spectrum of reference.
ref_spec = shaasp.cola_stft(ref_sig, nfft, hop, wind, single_sided);

% If ReTM is estimated as a shot-time correlation but is different size as
% the target recording, just revert back to an average time-invariant ReTM.
if size(retm_ab, 2) ~= size(ref_spec, 2)
    retm_ab = shaasp.spec_time_averaging(retm_ab, [], "all");
end

% Esimate noise spectrum at target by filtering reference with noise's ReTM.
est_spec = lastpagemtimes(retm_ab, ref_spec);  % [f,t,a]

% ISTFT to estimate noise signal at target microphones.
est_sig = shaasp.cola_istft(est_spec, nfft, hop);  % [samples, a]
est_sig(isnan(est_sig)) = 0;
est_sig = est_sig(1 : size(tgt_sig,1), :);

% Subtractive denoise target signal.
denoised_sig = tgt_sig - est_sig;

return

%% Subfunctions.

function [trim_sig] = trim_ref_to_tgt(ref_sig, tgt_sig)
    ref_nsamples = size(ref_sig, 1);
    tgt_nsamples = size(tgt_sig, 1);
    
    if (ref_nsamples < tgt_nsamples)
        trim_sig = [ref_sig; zeros(tgt_nsamples - ref_nsamples, 1)];
    elseif (ref_nsamples >= tgt_nsamples)
        trim_sig = ref_sig(1:tgt_nsamples, :);
    end
end

function [inv_x] = lastpagepinv(x)
    % Pagepinv applies to first 2 dimensions, but my mic channels are on
    % last 2 dimensions, so need permutes.
    n = ndims(x);
    x = permute(x, circshift((1:n), 2));  % [f,t,a,b] -> [a,b,f,t]
    % inv_x = pagepinv(x);  % If you have MATLAB 2024a.
    inv_x = shaasp.pagepinv(x);
    inv_x = permute(inv_x, circshift((1:n), -2));  % [b,a,f,t] -> [f,t,b,a]
end

function [z] = lastpagemtimes(x, y)
    n = ndims(x);
    x = permute(x, circshift((1:n), 2));  % [a,b,f,t]
    y = permute(y, circshift((1:n), 2));  % [b,c,f,t]
    z = pagemtimes(x, y);  % [a,c,f,t]
    z = permute(z, circshift((1:n), -2));  % [f,t,a,c]
end