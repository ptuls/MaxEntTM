function [cc,c,d] = clustcoeff(a)
%% CLUSTCOEFF - clustering coefficients of a symmetric adjacency matrix
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de

% test symmetry
issym = size(a,1)==size(a,2) && isempty(find(a-a'));
if ~issym
    error('CLUSTCOEFF requires symmetric input matrix');
end

% local clustering coefficient
d = node_degree(a); % degree of each node
cn = diag(a*triu(a)*a); % number of triangles for each node
c = zeros(size(d));
c(d>1) = 2*cn(d>1)./(d(d>1).*(d(d>1)-1));
%cc = sum(c.*d)/sum(d); % clustering coefficient 1
cc = mean(c(d>1));

if nargout==2
    % global clustering coefficient
    b = triu(a*a,1);
    c2 = sum(b); % number of connected triples
    c1 = sum(b.*a); % 3x number of triangles
    varargout{2} = full(sum(c1)/sum(c2)); % clustering coefficient 2
end