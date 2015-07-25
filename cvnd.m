function cv = cvnd(A)
%CVND Calculates coefficient of variation of the node degree (CVND)
%   
%   cv = cvnd(A) computes CVND given the adjacency matrix
%
%   [m,sd] = link_length(node_distances,adjacency) also returns the
%   standard deviation of the link lengths

% file:     cvnd.m, (c) Paul Tune, Fri Jan 09 2015
% created: 	Fri Jan 09 2015
% author:  	Paul Tune 
% email:   	paul.tune@adelaide.edu.au

degree = node_degree(A);
cv = std(degree)./mean(degree);