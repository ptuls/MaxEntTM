function x = peak_mean_cycle(freq,N,m,peakmean,troughmean)
%PEAK_MEAN_CYCLE Generate a simple sinusoid with specified frequency,
%   length, mean and peak-to-mean ratio
%
%   x = peak_mean_cycle(freq,N,m,peakmean) generates sinusoid via a 
%   sine function with frequency freq, length N, mean m and peak-to-mean
%   ratio peakmean. The generated signal has the form
%   
%   x = m*(peakmean-1)*sin(2*pi*freq*linspace(0,1,N))+m
%
%   Note that if the mean m is specified as 0, then 
%
%   x = peakmean*sin(2*pi*freq*linspace(0,1,N))
%
%   The peak-to-mean ratio must be greater than or equal to 1 (1 resulting
%   in a constant, flatlined signal).
%
%   x = peak_mean_cycle(freq,N,m,peakmean) generates a signal with  
%   a peak-to-mean ratio peakmean and a trough-to-mean ratio troughmean. 
%   This is because generated traffic may not be zero, which is problematic
%   if a pure sine function was used. Intervals of x below the mean now
%   have a drawdown determined by troughmean.
%
%   Note that the trough-to-mean ratio must be less than or equal to the
%   peak-to-mean ratio.
%

%
% file:     peak_mean_cycle.m, (c) Paul Tune, Thu Jun 04 2015
% created: 	Thu Jun 04 2015 
% author:  	Paul Tune 
% email:   	paul.tune@adelaide.edu.au

% modified: Fri Jun 12 2015
% Added trough-to-mean ratio to prevent zero traffic.

if peakmean < 1
    error('Peak-to-mean ratio must be greater than or equal to 1');
end

% what if the mean is zero?
if m == 0
    % take peak-to-mean ratio as the multiplicative factor
    x = peakmean*sin(2*pi*freq*linspace(0,N-1,N));
    return;
end

% know the peak is 1
x = m*(peakmean-1)*sin(2*pi*freq*linspace(0,N-1,N))+m;

if nargin > 4
    if troughmean > peakmean
        error('Trough-to-mean ratio must be less than or equal to peak-to-mean ratio');
    end
    if troughmean == 1
        return;     % no effect on trough
    end
    y = m*(troughmean-1)*sin(2*pi*freq*linspace(0,N-1,N))+m;
    x(x < m) = y(y < m);
end