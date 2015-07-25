function TM = ncim(row,col,total,sigmasq)
%NCIM Generate a series of traffic matrices via the Non-stationary 
% Conditionally Independent Model (NCIM). 
%
% TM = ncim(row,column,total,sigmasq) returns the series of traffic 
% matrices TM based on 
%
%  1. row: average total outgoing traffic from PoPs in each time interval
%  2. col: average total incoming traffic into PoPs in each time interval
%  3. total: total_traffic in each time interval 
%  4. sigmasq: noise variance \sigma^2 for input into the truncated normal
%     distribution
%
%  The script generates M traffic matrices, where M is the length of
%  total, with the spatial size determined by the number of rows 
%  row and number of columns of col. Note that row and col must have their
%  lengths equal to that of total, i.e. they must be defined on M time
%  intervals.
%
%  Requires the HMC_exact script by Ari Pakman. Download from 
%  https://github.com/aripakman/hmc-tmg

%
% file:     ncim.m, (c) Paul Tune, Thu May 28 2015
% created: 	Thu May 28 2015 
% author:  	Paul Tune 
% email:   	paul.tune@adelaide.edu.au

% sanity checks
if isempty(total)
    error('No total traffic provided');
end

[s,t] = size(total);
if s > 1 && t > 1
    error('Total traffic should not be a matrix');
end

if s > 1
    total = total';
end

if min(min(row)) < 0
    error('Row entries must be non-negative');
end

if min(min(col)) < 0
    error('Column entries must be non-negative');
end

[Nrow,n_TMs] = size(row);
[Ncol,m_TMs] = size(col);
if Nrow ~= Ncol
    error('Column length mismatch to row length');
end

if n_TMs ~= m_TMs
    error('Column temporal length mismatch to row temporal length');
end

if length(total) ~= n_TMs
    error('Total traffic length mismatch to row/column temporal length');
end

if sigmasq < 0
    error('Noise variance must be non-negative');
end

if min(modulated_total) < 0
    error('Total traffic must be non-negative');
end

%% construct NCIM
TM = zeros(Nrow,Ncol,n_TMs);
for i=1:n_TMs 
    % fanouts
    pU = row(:,i)/total(i);
    pV = col(:,i)/total(i);
    
    % construct gravity-like model in each time slice
    U = HMC_exact(eye(Nrow), zeros(Nrow,1), eye(Nrow)*sigmasq, pU, true, 2, pU);
    V = HMC_exact(eye(Ncol), zeros(Ncol,1), eye(Ncol)*sigmasq, pV, true, 2, pV);

    TM(:,:,i) =  total(i)*U(:,2)*V(:,2)';
end