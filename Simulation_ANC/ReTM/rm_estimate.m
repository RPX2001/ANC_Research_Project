function [rm, rmi] = rm_estimate(tgt_sig, ref_sig, options)
% RM_ESTIMATE - Estimate RM mapping: TARGET = RM * REFERENCE (frequency-wise LS)
% Outputs:
%   rm  - [f,1,A,B] frequency-domain RM (A target channels, B reference channels)
%   rmi - [sample,1,A,B] time-domain impulse responses via IFFT

arguments
    tgt_sig (:,:)
    ref_sig (:,:)
    options.fs = []
    options.wlen = []
    options.nfft = []
    options.hop = []
    options.reg = 1e-6
end

fs   = options.fs;
wlen = options.wlen;
nfft = options.nfft;
hop  = options.hop;
reg  = options.reg;

if isempty(wlen)
    if ~isempty(fs)
        wlen = fs / 1;
    else
        wlen = 4096;
    end
end

if isempty(nfft)
    nfft = 2^nextpow2(2 * wlen);
end

if isempty(hop)
    hop = wlen / 4;
end

if size(tgt_sig, 1) ~= size(ref_sig, 1)
    error('Expect input signals to be the same length.');
end

ola_method = 'ola';
wind = shaasp.cola_window(wlen, hop, ola_method);
single_sided = true;

% STFT
E = shaasp.cola_stft(tgt_sig, nfft, hop, wind, single_sided); % [F,T,A]
M = shaasp.cola_stft(ref_sig, nfft, hop, wind, single_sided); % [F,T,B]

[Fbins, Tframes, A] = size(E);
B = size(M,3);

rm = zeros(Fbins, 1, A, B);

% Frequency-wise multi-channel ridge least squares:
% minimize ||E_f - M_f * W_f^T||_F^2
% W_f = (M_f^H M_f + reg*I)^(-1) M_f^H E_f   then transpose to [A x B]
for f = 1:Fbins
    Mf = squeeze(M(f,:,:));  % [T x B]
    Ef = squeeze(E(f,:,:));  % [T x A]

    if isvector(Mf), Mf = reshape(Mf, [], B); end
    if isvector(Ef), Ef = reshape(Ef, [], A); end

    G = (Mf' * Mf) + reg * eye(B);
    W = G \ (Mf' * Ef);      % [B x A]

    rm(f,1,:,:) = permute(W.', [3 4 1 2]); % [1,1,A,B]
end

rmi = ifft(rm, nfft, 1, 'symmetric');
end
