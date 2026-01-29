function [retf, reir] = retf_estimate(tgt_sig, ref_sig, options)
% RETF_ESTIMATE - Estimate the ReTF between two signals.
% Estimates the ReTF for the mapping: TARGET = R * REFERENCE.
%
% Syntax:  [retf, reir] = retf_estimate(tgt_sig, ref_sig, options);
%          [retf, reir] = retf_estimate(tgt_sig, ...
%                                       ref_sig, ...
%                                       'fs', fs, ...
%                                       'wlen', wlen, ...
%                                       'nfft', nfft, ...
%                                       'hop', hop);
%
% Inputs:
%   tgt_sig - [samples, chn=1] Target signal (ReTF maps ref -> target).
%   ref_sig - [samples, chn=1] Reference signal.
%   - Optional name.value inputs - 
%   fs - Sampling frequency of both tgt_sig and ref_sig.
%   wlen - stft window length.
%   nfft - stft nfft (with zero padding) size.
%   hop - stft hop size.
%
% Outputs:
%   retf - Estimated ReTF (frequency domain).
%   reir - IFFT(retf) with correct nfft size.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: example_retf_denoise
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 10-Oct-2024
% Last revision: 10-Oct-2024

arguments
    tgt_sig (:,1)
    ref_sig (:,1)
    options.fs = []
    options.wlen = []
    options.nfft = []
    options.hop = []
end

fs = options.fs;
wlen = options.wlen;
nfft = options.nfft;
hop = options.hop;

if isempty(wlen) && ~isempty(fs)
    wlen = fs / 1;  % Default 1 second window.
else
    wlen = 4096;
end

if isempty(nfft)
    nfft = 2^nextpow2(2 * wlen);  % Default nfft is 2*Wlen.
end

if isempty(hop)
    hop = wlen / 4;
end

if length(tgt_sig) ~= length(ref_sig)
    error('Expect input signals to be the same length.');
end

ola_method = 'ola';
wind = shaasp.cola_window(wlen, hop, ola_method);

% STFT signals.
tgt_spec = shaasp.cola_stft(tgt_sig, nfft, hop, wind);
ref_spec = shaasp.cola_stft(ref_sig, nfft, hop, wind);

% CPSD spectrums.
cpsd_12 = shaasp.spec_to_crossspec(tgt_spec, ref_spec);
cpsd_22 = shaasp.spec_to_crossspec(ref_spec);

% Time averaging to get spatial correlation.
corr_12 = mean(cpsd_12, 2); % shaasp.spec_time_averaging(cpsd_12, [], "all");
corr_22 = mean(cpsd_22, 2); % shaasp.spec_time_averaging(cpsd_22, [], "all");

% Estimate ReTF.
retf = corr_12 ./ corr_22;
reir = real(ifft(retf, nfft, 1, "symmetric"));

end

