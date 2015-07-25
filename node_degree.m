function degree = node_degree(A)
%NODE_DEGREE Computes the node degree of an undirected graph
%   
%   degree = node_degree(A) computes node degree of an undirected graph,
%   given adjacency matrix A
%

% file:     node_degree.m, (c) Paul Tune, Fri Jan 09 2015
% created: 	Fri Jan 09 2015
% author:  	Paul Tune 
% email:   	paul.tune@adelaide.edu.au

% check if A is symmetric
if A == (A+A')/2
    error('Adjacency matrix is not for an undirected graph');
end

degree = sum(A,2);