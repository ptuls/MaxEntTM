%%GAmaxent_minimax
% Testing COLD on a series of traffic matrices via the minimax scheme. 
% Locations are generated once, loaded from file example_topo.mat. 
% Modulated gravity model was used
% and Fourier coefficients extracted from Abilene's total traffic.
% Abilene's traffic in 5 minute intervals is found in tm_real.mat. 28
% Fourier coefficients were chosen via nMDL, computed beforehand.
%
%
% file:     GAmaxent_minimax.m, (c) Paul Tune, Fri Jan 09 2015
% created: 	Fri Jan 09 2015
% author:  	Paul Tune 
% email:   	paul.tune@adelaide.edu.au
% 


load tm_real
load example_topo


%%%%%%%%%%%%%%%%%%%%% PARAMETER SETTINGS FOR EXPERIMENT AND COLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trial number
trial = 1;

%% Noise strength
low = 1;

%% Parameters for COLD
num_nodes = 12;
width = 10;
height = 10;
node_population_param = 20;
pop_distn = 'exp';

% optimisation parameters
parameters = struct();
parameters.k0 = 5;
parameters.k1 = 1;
parameters.k2 = k2;
parameters.k3 = k3;

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
preload_locations = 1;

% generate new locations if not preloading locations
if (~preload_locations)
    placesX = unifrnd(0,width,[num_nodes 1]);   % conditional 2D Poisson
    placesY = unifrnd(0,height,[num_nodes 1]);
    places_pop = random(pop_distn,node_population_param,num_nodes,1); % exponential population model
    node_map = [placesX placesY places_pop];
    node_distances = zeros(num_nodes);
end
parameters.node_distances = node_distances;

%%%%%%%%%%%%%%%%%%%%% MAXENT MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modulated gravity model
n_TMs = 2016;
tm_tensor = reshape(tm_real,N,N,n_TMs);
r = squeeze(sum(tm_tensor,2));
c = squeeze(sum(tm_tensor,1));

% total traffic in each interval
T = sum(tm_real,1);

% mean row, column and total sums
mean_r = mean(r,2);
mean_c = mean(c,2);
mean_T = mean(T);

% construct k-sparse version of total traffic
% extract Fourier coefficients and form k-sparse version
fT = fft(T-mean_T);
[fTmag,indf] = sort(abs(fT(1:(n_TMs/2))),2,'descend');

% choose k coefficients based on nMDL
k = 23;
fs = zeros(1,n_TMs/2);
fs(indf(1:k)) = fT(indf(1:k));
fsparse = [fs fT(n_TMs/2+1) conj(fs(end:-1:2))];
norm_T = real(ifft(fsparse));

%% modulated gravity based on sparse Fourier, nMDL selected coefficients
tic
X = zeros(num_nodes,num_nodes,n_TMs);
Tk = mean(T);
rk = mean(r,2);
ck = mean(c,2);

if low
	sigmasq = 1e-5;
else
	sigmasq = 10;
end
pU = rk/Tk;
pV = ck/Tk;
U = HMC_exact(eye(N), zeros(N,1), eye(N)*sigmasq/Tk^2, pU, true, 2, pU);
V = HMC_exact(eye(N), zeros(N,1), eye(N)*sigmasq/Tk^2, pV, true, 2, pV);
M0 = (Tk+norm_T)/Tk;
M = HMC_exact(eye(n_TMs),zeros(n_TMs,1),sigmasq/Tk^2,M0',true,2,M0');

G = Tk*U(:,2)*V(:,2)';
X = reshape(reshape(G,N^2,1)*M(:,2)',N,N,n_TMs);
toc

%% Minimax
% minimum: 1608
% maximum: 308
M = n_TMs;
% M = 288; % one day
% M= 10;
parameters.demand = X(:,:,1:M);
[topo,best_cost] = coldGA(parameters);

[gcc,gc] = clustcoeff(topo.adjacency);

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
% minimax
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

if low
	filename = 'minimax_';
else
	filename = 'minimax_s10_';
end
params = ['k2_' num2str(parameters.k2,'%.0e') '_k3_' ...
    num2str(parameters.k3,'%.0e') '_trial_' num2str(trial,'%d') '.mat'];

filename = [filename params];

while exist(filename)
	trial = trial+1;
if low
	filename = 'minimax_';
else
	filename = 'minimax_s10_';
end
params = ['k2_' num2str(parameters.k2,'%.0e') '_k3_' ...
    num2str(parameters.k3,'%.0e') '_trial_' num2str(trial,'%d') '.mat'];
filename = [filename params];

	if (trial > 20 )
		break;
	end
end

save(filename,'results','topo');
fprintf('Saving file %s\n',filename);
