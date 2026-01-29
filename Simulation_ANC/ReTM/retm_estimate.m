function [retm, reim] = retm_estimate(tgt_sig, ref_sig, options)
% RETM_ESTIMATE - Estimate the ReTM between two multi-channel signals.
% Estimates the ReTM for the mapping: TARGET = R * REFERENCE. 
%
% Syntax:  [retm, reim] = retm_estimate(tgt_sig, ref_sig, options);
%          [retm, reim] = retm_estimate(tgt_sig, ...
%                                       ref_sig, ...
%                                       'fs', fs, ...
%                                       'wlen', wlen, ...
%                                       'nfft', nfft, ...
%                                       'hop', hop);
%
% Inputs:
%   tgt_sig - [samples, chn_A] Target signals (ReTM maps ref -> target).
%   ref_sig - [samples, chn_B] Reference signals.
%   - Optional name.value inputs - 
%   fs - Sampling frequency of both tgt_sig and ref_sig.
%   wlen - stft window length.
%   nfft - stft nfft (with zero padding) size.
%   hop - stft hop size.
%
% Outputs:
%   retm - [f,1,A,B] Estimated ReTM (frequency domain).
%   reim - [sample,1,A,B] IFFT(retm) with correct nfft size.
%
% Other m-files required: none
% Subfunctions: lastpagepinv(), lastpagemtimes().
% MAT-files required: none
%
% See also: example_retm_denoise
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 10-Oct-2024
% Last revision: 10-Oct-2024

arguments
    tgt_sig (:,:)
    ref_sig (:,:)
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

if size(tgt_sig, 1) ~= size(ref_sig, 1)
    error('Expect input signals to be the same length.');
end

ola_method = 'ola';
wind = shaasp.cola_window(wlen, hop, ola_method);
single_sided = true;  % Only do half-spectrum for better performance.

% STFT signals.
tgt_spec = shaasp.cola_stft(tgt_sig, nfft, hop, wind, single_sided);
ref_spec = shaasp.cola_stft(ref_sig, nfft, hop, wind, single_sided);

% CPSD spectrums.
cpsd_aa = shaasp.spec_to_crossspec(tgt_spec);  % [f,t,a,a]
cpsd_ba = shaasp.spec_to_crossspec(ref_spec, tgt_spec);  % [f,t,b,a]

% Time averaging to get spatial correlation matrix.
corr_aa = mean(cpsd_aa, 2); % shaasp.spec_time_averaging(cpsd_aa, [], "all");
corr_ba = mean(cpsd_ba, 2); % shaasp.spec_time_averaging(cpsd_ba, [], "all");

% Estimate inverse of corr_ba. 
corr_ba_inv = lastpagepinv(corr_ba);  % [f,1,a,b]

% Estimate ReTM = corr_aa * (corr_ba)^-1
retm = lastpagemtimes(corr_aa, corr_ba_inv);  % [f,1,a,b]
reim = ifft(retm, nfft, 1, 'symmetric');  % [sample,1,a,b]

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