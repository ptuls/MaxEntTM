function [TM,G] = modulated_gravity(mean_row,mean_col,modulated_total,sigmasq,sigmasqt)
%MODULATED_GRAVITY Generate a series of traffic matrices via the
%  Modulated Gravity Model (MGM). This model assumes that the average
%  spatial traffic matrix and its variance is known.
%
%  TM = modulated_gravity(mean_row,mean_col,modulated_total,sigmasq,sigmasqt)
%  returns a series of traffic matrices TM from the modulated gravity model
%  given:
%
%  1. mean_row: average total outgoing traffic from PoPs
%  2. mean_col: average total incoming traffic to PoPs
%  3. modulated_total: total_traffic over the time interval 
%  4. sigmasq: noise variance \sigma^2 for input into the truncated normal
%     distribution
%  5. sigmasqt: noise variance for the temporal component, i.e. variance of
%     total traffic (optional). If this is not set, the default is
%     sigmasqt = sigmasq
%
%  The script generates M traffic matrices, where M is the length of
%  modulated_total, with the spatial size determined by mean_row and
%  mean_col.
%
%  [TM,G] = modulated_gravity(mean_row,mean_col,modulated_total,sigmasq,sigmasqt)
%  also returns the gravity model G for reference.

%  Requires the HMC_exact script by Ari Pakman. Download from 
%  https://github.com/aripakman/hmc-tmg

%
% file:     modulated_gravity.m, (c) Paul Tune, Wed May 27 2015
% created: 	Wed May 27 2015 
% author:  	Paul Tune 
% email:   	paul.tune@adelaide.edu.au

% sanity checks
if nargin < 5
    sigmasqt = sigmasq;
end

if isempty(modulated_total)
    error('No modulated total traffic provided');
end

[s,t] = size(modulated_total);
if s > 1
    modulated_total = modulated_total';
end

if s > 1 && t > 1
    error('Modulated signal should not be a matrix');
end

if min(mean_row) < 0
    error('Row means must be non-negative');
end

if min(mean_col) < 0
    error('Column means must be non-negative');
end

N = length(mean_row);
if length(mean_col) ~= N
    error('Column means length mismatch to row means length');
end

if sigmasq < 0
    error('Noise variance must be non-negative');
end

% force row and column means to the correct configuration
[~,t] = size(mean_row);
if t > 1
    mean_row = mean_row';
end

[~,t] = size(mean_col);
if t > 1
    mean_col = mean_col';
end

if min(modulated_total) < 0
    error('Total traffic must be non-negative');
end
% setup necessary auxiliary parameters: number of TMs to generate and the
% mean of the total traffic from the modulated total traffic signal
n_TMs = length(modulated_total);
mean_total = mean(modulated_total);

% average fanout
pU = mean_row/mean_total;
pV = mean_col/mean_total;

% Synthesize
% generate truncated normal random variables: we generate two samples
% each time because HMC_exact needs to burn-in 
U = HMC_exact(eye(N), zeros(N,1), eye(N)*sigmasq/mean_total^2, pU, true, 2, pU);
V = HMC_exact(eye(N), zeros(N,1), eye(N)*sigmasq/mean_total^2, pV, true, 2, pV);

% modulate mean with modulated reference total
normalized_mean = modulated_total/mean_total;
modulated = HMC_exact(eye(n_TMs),zeros(n_TMs,1),sigmasqt/mean_total^2,normalized_mean',true,2,normalized_mean');
disp(modulated);
% gravity model
G = mean_total*U(:,2)*V(:,2)';

% construct modulated gravity model
% in 3D array form
TM = reshape(reshape(G,N^2,1)*modulated(:,2)',N,N,n_TMs);