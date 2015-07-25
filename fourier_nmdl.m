function [Snmdl,Knmdl] = fourier_nmdl(S)
%FOURIER_NMDL Computing a frequency sparse signal from original signal based on 
% the normalized Minimum Description Length (nMDL) information
% criterion.
%
% Likelihood is based on multivariate normal distribution, uncorrelated
% noise, once the signal has been centered. A direct likelihood
% using a multivariate truncated normal distribution is also possible, but
% no closed form may exist for the nMDL criterion.
%
% Criterion used here is based on the approximation derived in Hansen and Yu, 
% "Model Selection and the Principle of Minimum Description Length", 
% Journal of the American Statistical Association, Vol. 96, No. 454 
% (June 2001), pp. 746-774.
%
% Snmdl = fourier_nmdl(S) returns sparse version of signal S based on 
% selecting Fourier coefficients via nMDL.
%
% [Snmdl,Knmdl] = fourier_nmdl(S) also returns Knmdl, the number of Fourier
% coefficients selected.

% file:     fourier_nmdl.m, (c) Paul Tune, Thu May 28 2015
% created: 	Thu May 28 2015 
% author:  	Paul Tune 
% email:   	paul.tune@adelaide.edu.au

[m,n] = size(S);
if m > 1 && n > 1
    error('Not a 1-dimensional signal');
end

if m > 1
    S = S'; % turn into row vector
end

if mod(length(S),2) ~= 0
    % means signal has odd length, add DC component
    S = [S mean(S)];
end
Slen = length(S);

% preprocessing: center signal by removing mean
Smean = mean(S);
fS = fft(S-mean(S));
[~,indf] = sort(abs(fS(1:(Slen/2))),2,'descend');

lim = Slen/4;
nMDL = zeros(1,lim);

% search for least number of Fourier coefficients based on nMDL
% unfortunately, brute force based at the moment
for k=1:lim
    fs = zeros(1,Slen/2);
    
    fs(indf(1:k)) = fS(indf(1:k));
    fsparse = [fs fS(Slen/2+1) conj(fs(end:-1:2))];
    
    x(k,:) = real(ifft(fsparse)) + Smean;
    
    % full Fourier, no time sparse component
    % each Fourier component has 2 DOF: magnitude and phase
    residue = norm(S-x(k,:));
    
    % number of parameters
    p = 2*k+1;
    
    % compute nMDL criterion: see Hansen and Yu, "Model Selection and the
    % Principle of Minimum Description Length", Journal of the American 
    % Statistical Association, Vol. 96, No. 454 (June 2001), pp. 746-774 
    R = residue^2/(Slen-p);
    nMDL(k) = 0.5*Slen*log(R) + 0.5*p*log((norm(S).^2 ...
        - residue^2)/(p*R)) + 0.5*log(Slen-p) - 1.5*log(p);
end

% extract Fourier coefficients and form sparse version from nMDL
% choose k coefficients based on nMDL
[~,Knmdl] = min(nMDL');

% reconstruct sparse signal
fs = zeros(1,Slen/2);
fs(indf(1:Knmdl(1))) = fS(indf(1:Knmdl(1)));
fsparse = [fs fS(Slen/2+1) conj(fs(end:-1:2))];
Snmdl = real(ifft(fsparse)) + Smean; % clean from machine noise