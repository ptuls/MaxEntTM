function [xs,fsparse] = freq_sparsify(x,K)
%FREQ_SPARSIFY Computing a frequency sparse signal with K coefficients from
% original signal.
%
% xs = freq_sparsify(x,K) returns the frequency sparse component
% xs with exactly K Fourier coefficients.
%
% [xs,fsparse] = freq_sparsify(x,K) also returns fsparse, the discrete
% Fourier transform of xs.

% file:     freq_sparsify.m, (c) Paul Tune, Thu May 28 2015
% created: 	Thu May 28 2015 
% author:  	Paul Tune 
% email:   	paul.tune@adelaide.edu.au


lx = length(x);
fx = fft(x-mean(x));
[fxmag,indf] = sort(abs(fx(1:(lx/2))),2,'descend');

fs = zeros(1,lx/2);
fs(indf(1:K)) = fx(indf(1:K));
fsparse = [fs fx(lx/2+1) conj(fs(end:-1:2))];
fsparse(1) = sum(x);

xs = real(ifft(fsparse));
fsmag = fxmag(1:K);