function [m,sd] = link_length(node_distances,adjacency)
%LINK_LENGTH Calculates the average and standard deviation of the 
%   link lengths
%   
%   m = link_length(node_distances,adjacency) computes the average link
%   length, given the node distances and adjacency matrix
%
%   [m,sd] = link_length(node_distances,adjacency) also returns the
%   standard deviation of the link lengths

% file:     link_length.m, (c) Paul Tune, Fri Jan 09 2015
% created: 	Fri Jan 09 2015
% author:  	Paul Tune 
% email:   	paul.tune@adelaide.edu.au

if size(node_distances) ~= size(adjacency)
    error('Node distance matrix and adjacency matrix not the same size');
end

% ensure non-zero entries set to 1
adjacency(abs(adjacency) > 0) = 1;

L = node_distances.*triu(adjacency);

m = mean(L(L > 0));
sd = std(L(L > 0));