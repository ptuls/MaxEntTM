%% CASE_STUDY_TOPOLOGY
% Testing COLD on a series of traffic matrices via the minimax scheme. 
% Locations are generated once, loaded from file example_topo.mat. 
% Modulated gravity model (MGM) was used
% and Fourier coefficients extracted from Abilene's total traffic via the 
% normalized Minimum Description Length (nMDL) information criterion.
% Abilene's traffic in 5 minute intervals is found in tm_real.mat. 28
% Fourier coefficients were chosen via nMDL, computed beforehand.
%
% Note: this generates the results for a selected instance of COLD's cost
% k2 and k3. All results are saved in a struct called results. The
% following are saved:
%
% 1. total_links: total number of links in the designed topology
% 2. average_degree: the average degree of PoPs in the designed topology
% 3. cvnd: coefficient of variations of the node degree
% 4. diameter: diameter of the topology
% 5. gcc: the global clustering coefficient
% 6. best_cost: the cost of the optimal topology found by COLD
% 7. mean_real: average link load of the topology if Abilene's is used
% 8. capacity_real: link capacity of topology using Abilene's traffic
% 9. mean_synthetic: average link load of the topology if Abilene's is used
% 10. capacity_synthetic: link capacity of topology using MGM's traffic
%
% These data can be used to do a variety of things. For instance, we can
% compare the topology's link capacity on real traffic versus synthetic
% trafic. We can also compute other quantities such as maximum utilization.

% file:     case_study_topology.m, (c) Paul Tune, Fri Jan 09 2015
% created: 	Fri Jan 09 2015
% author:  	Paul Tune 
% email:   	paul.tune@adelaide.edu.au
% 

clear;
% load Abilene traffic matrices and network topology
load tm_real

%% Parameters for traffic matrix generation
low = 1; % here we toggle between \sigma^2 = 10^-5 and 10

% number of traffic matrices to generate
% minimum: 1608
% maximum: 308
M = n_TMs;
% M = 288; % one day

% set to 0 if new topology is to be generated and the example topology is
% not to be used
preload_locations = 1;

%% Parameters for COLD
num_nodes = 12;
width = 10;
height = 10;
node_population_param = 20;
pop_distn = 'exp';

if preload_locations
    load example_topo;
    % overwrite existing parameters
    num_nodes = 12;
    width = 10;
    height = 10;
    node_population_param = 20;
    pop_distn = 'exp';
end

% optimisation parameters
parameters = struct();
parameters.k0 = 5;
parameters.k1 = 1;
parameters.k2 = 1;
parameters.k3 = 0.5;

% internal GA parameters
parameters.num_nodes = num_nodes;
parameters.width = width;
parameters.height = height;
parameters.node_population_param = node_population_param;
parameters.pop_distn = pop_distn;
parameters.num_chromosomes=100;
parameters.num_generations=100;
parameters.crossovera=2;
parameters.crossoverb=10;
parameters.mutate_prob = 0.02;
parameters.start_links = 2;
parameters.mutate_fn = @()geornd([0.4,0.4]);
parameters.num_saved_chromosomes = 50;
parameters.fast_mode =1;
parameters.hub_mut_rate = 0.05;
parameters.crossover_rate = 0.3;
parameters.criterion = 'minimax';   % minimax scheme

%% Parameters for PoP locations

% generate new locations if not preloading locations
if (~preload_locations)
    placesX = unifrnd(0,width,[num_nodes 1]);   % conditional 2D Poisson
    placesY = unifrnd(0,height,[num_nodes 1]);
    places_pop = random(pop_distn,node_population_param,num_nodes,1); % exponential population model
    node_map = [placesX placesY places_pop];
    node_distances = zeros(num_nodes);
end
parameters.node_distances = node_distances;

%%%%%%%%%%%%%%%%%%%%% MAXENT MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modulated gravity model
n_TMs = 2016;
tm_tensor = reshape(tm_real,N,N,n_TMs);
r = squeeze(sum(tm_tensor,2));
c = squeeze(sum(tm_tensor,1));

% total traffic in each time interval
T = sum(tm_real,1);

% compute mean
mean_row = mean(r,2);
mean_col = mean(c,2);
mean_T = mean(T);

% extract Fourier coefficients and form sparse version, where Fourier
% coefficients are selected via the nMDL criterion
[Tsparse,~,Ksparse] = fourier_time_nmdl(T);

fprintf('Total traffic - Fourier coefficients selected: %d\n', Ksparse);

%% modulated gravity based on sparse Fourier, nMDL selected coefficients
mean_Tsparse = mean(Tsparse);
if low
	sigmasq = 1e-5;
else
	sigmasq = 10;
end

tic
TMmod_gravity = modulated_gravity(mean_row,mean_col,Tsparse,sigmasq);
toc

%% Minimax

parameters.demand = X(:,:,1:M);
[topo,best_cost] = coldGA(parameters);

[gcc,gc] = clustcoeff(topo.adjacency);

% store results in struct
results = struct();
results.total_links = sum(sum(triu(topo.adjacency)));
results.average_degree = mean(node_degree(topo.adjacency));
results.cvnd = cvnd(topo.adjacency);
results.diameter = diameter(topo.adjacency);
results.gcc = gcc;
results.best_cost = best_cost;

[m,sd] = link_length(node_distances,topo.adjacency);
results.mean_link_length = m;
results.std_link_length = sd;

%% Calculate maximum load of links in each topology over the
%  M TMs to see which requires less capacity

% follow the order as done in COLD
[distance_order,ix] = sort(reshape(node_distances,[],1));
distance_order=zeros(size(ix));
for ind = 1:length(ix)
    distance_order(ix(ind))=ind;
end
distance_order = reshape(distance_order,num_nodes,num_nodes);

BIG = num_nodes*(max(max(node_distances))+1);
dists1 = zeros(num_nodes);

%% Compute path distances for all topologies using Floyd-Warshall

dists1(topo.adjacency>0) = node_distances(topo.adjacency>0);
dists1(topo.adjacency==0) = BIG*distance_order(topo.adjacency==0);
% do the routing
[path_distances1, pred1] = floyd_apsp(dists1);

% assume we find best shortest path routing per measurement interval
% we calculate the link loads for each topology based on best routing
ospf = 1; % ospf = 0 ==> disable equal-cost path splitting
capacity_real = zeros(num_nodes,num_nodes);
capacity_synthetic = zeros(num_nodes,num_nodes);

load1 = zeros(num_nodes, num_nodes, n_TMs);
load2 = zeros(num_nodes, num_nodes, n_TMs);
for i=1:n_TMs
    load1(:,:,i) = calc_link_load(path_distances1, dists1, pred1, tm_tensor(:,:,i), ospf);
    load2(:,:,i) = calc_link_load(path_distances1, dists1, pred1, X(:,:,i), ospf);
    
    capacity_real = max(capacity_real,load1(:,:,i));
    capacity_synthetic = max(capacity_synthetic,load2(:,:,i));
end

results.mean_real = mean(load1(load1 > 0));
results.capacity_real = capacity_real;
results.mean_synthetic = mean(load2(load2 > 0));
results.capacity_synthetic = capacity_synthetic;

%% Plot topology
fontsize = 20;
axeswidth = 3;
markersize = 20;
thick = 2;

figure(1)
gplot(topo.adjacency,node_map(:,1:2));
hold on
plot(node_map(:,1),node_map(:,2),'.','color','r', 'markersize',markersize);
hold off
set(findall(gcf,'type','line'),'LineWidth',thick)

hx = xlabel('x');
hy = ylabel('y','Rotation',0);
axes_handle = gca;
set(axes_handle, 'FontSize', fontsize,'LineWidth',axeswidth);
set(hx, 'FontSize', fontsize);
set(hy, 'FontSize', fontsize);