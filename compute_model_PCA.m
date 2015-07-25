%% COMPUTE_MODEL_PCA Compute the principle component analysis of the
%  Modulated Gravity Model (MGM) and the Non-stationary Conditionally
%  Independent Model (NCIM). 

% file:     maxent_pca.m, (c) Paul Tune, Mon Dec 22 2014
% created: 	Mon Dec 22 2014 
% author:  	Paul Tune 
% email:   	paul.tune@adelaide.edu.au

clear;

%% noise parameter \sigma^2 for both models
sigmasq = 1e-2;

%% get Abilene data
% load Abilene traffic matrices; we will extract the Fourier coefficients
% from these series of traffic matrices
load tm_real

% Note: tm_real already contains the spatial size (N) and temporal size
% (n_TMs) parameters

% turn into 3D array and compute row and column marginals
tm_tensor = reshape(tm_real,N,N,n_TMs);
r = squeeze(sum(tm_tensor,2));  % row marginals
c = squeeze(sum(tm_tensor,1));  % column marginals

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
tic
TMmod_gravity = modulated_gravity(mean_row,mean_col,Tsparse,sigmasq);
toc

%% non-stationary conditionally independent model (NCIM) with delayed rows,
%  columns to represent time zone shifts

% need to construct the row and column spatiotemporal signals
% choose New York's as reference
idNY = 9;
rs = zeros(N,n_TMs);
cs = zeros(N,n_TMs);
rs(idNY,:) = freq_sparsify(r(idNY,:),Ksparse); % nMDL calculated optimal choice of frequencies
cs(idNY,:) = freq_sparsify(c(idNY,:),Ksparse);

% each measurement interval is 5 minutes, so need to convert hours to
% minutes
m = 60/5;

% delay due to time zones (everywhere except NY)
D = m*[0 0 1 2 1 0 1 3 3 3 0];
location = [1:8 10:N]; % everywhere except NY

% mean_rs = mean(rs,2);
% mean_cs = mean(cs,2);

% build constraints: note we need constraints for every time interval
rscale = mean_row/mean_row(idNY);
cscale = mean_col/mean_col(idNY);
for i=1:length(location)
    rs(location(i),:) = circshift(rscale(location(i))*rs(idNY,:),[0 D(i)]);
    cs(location(i),:) = circshift(cscale(location(i))*cs(idNY,:),[0 D(i)]);
end
Ts = sum(rs,1);

tic
TMncim = ncim(rs,cs,Ts,sigmasq);
toc

%% Compute traffic decomposition via PCA
x_real = reshape(tm_real,N^2,n_TMs);
x_mod_gravity = reshape(TMmod_gravity,N^2,n_TMs);
x_ncim = reshape(TMncim,N^2,n_TMs);

mean_x_real = mean(x_real,2);
mean_x_mod_gravity = mean(x_mod_gravity,2);
mean_x_ncim = mean(x_ncim,2);

coeff_real = princomp(x_real);
coeff_mod_gravity = princomp(x_mod_gravity);
coeff_ncim = princomp(x_ncim);

spatial_real = svd(x_real);
[spatial_mod_gravity] = eig(x_mod_gravity*x_mod_gravity'-mean_x_mod_gravity *mean_x_mod_gravity');
spatial_ncim = eig(x_ncim*x_ncim'-mean_x_ncim*mean_x_ncim');

temporal_real = eig(x_real'*x_real);
temporal_mod_gravity = eig(x_mod_gravity'*x_mod_gravity - mean_x_mod_gravity'*mean_x_mod_gravity);
temporal_ncim = eig(x_ncim'*x_ncim - mean_x_ncim'*mean_x_ncim);


% remove machine noise
spatial_real(spatial_real < 0) = 0;
spatial_mod_gravity(spatial_mod_gravity < 0) = 0;
spatial_ncim(spatial_ncim < 0) = 0;

% remove machine noise
temporal_real(temporal_real < 0) = 0;
temporal_mod_gravity(temporal_mod_gravity < 0) = 0;
temporal_ncim(temporal_ncim < 0) = 0;

%% Parameters for plotting
dgreen = [0.05 0.5 0.05];
dgrey = [0.65 0.65 0.65];
fontsize = 18;
print_figures = 0;  % set to print figures into .eps, or not
paper_position = [0 0 14 12];
device = '-depsc';
thick_line = 1.1;
line = 1;
outfile = 'maxent_pca_';

%% Plot total sum
figure(1)
plot(sum(tm_real,1),'color',dgrey,'linewidth',thick_line);
hold on
plot(squeeze(sum(sum(TMmod_gravity))),'color',dgrey);
plot(squeeze(sum(sum(TMncim))),'color',dgreen);
hold off
set(gca,'fontsize', fontsize);
xlim([0,n_TMs]);
xlabel('t_k');
ylabel('Gbps');

%% Plot first three coeffs of real, modulated gravity model and 
%  non-stationary conditionally independent model
figure(10)
plot(coeff_real(:,1),'color',dgrey,'linewidth',thick_line);
hold on
plot(coeff_mod_gravity(:,1),'color',dgrey);
plot(coeff_ncim(:,1),'color',dgreen);
hold off
set(gca,'fontsize', fontsize);
xlim([0,n_TMs]);
xlabel('t_k');
ylabel('Coefficient magnitude');

figure(11)
plot(coeff_real(:,2),'color',dgrey,'linewidth',thick_line);
hold on
plot(coeff_mod_gravity(:,2),'color',dgrey);
plot(coeff_ncim(:,2),'color',dgreen);
hold off
set(gca,'fontsize', fontsize);
xlim([0,n_TMs]);
xlabel('t_k');
ylabel('Coefficient magnitude');

figure(12)
plot(coeff_real(:,3),'color',dgrey,'linewidth',thick_line);
hold on
plot(coeff_mod_gravity(:,3),'color',dgrey);
plot(coeff_ncim(:,3),'color',dgreen);
hold off
set(gca,'fontsize', fontsize);
xlim([0,n_TMs]);
xlabel('t_k');
ylabel('Coefficient magnitude');

%% Compare coeffs between real TM and models
figure(13)
id = 5;
L1 = length(coeff_real(:,id));
L2 = length(coeff_mod_gravity(:,id));
L3 = length(coeff_ncim(:,id));

semilogy(sort(coeff_real(:,id)),(1:L1)/L1,'color',dgrey,'linewidth',thick_line);
hold on
semilogy(sort(coeff_mod_gravity(:,id)),(1:L2)/L2,'r');
semilogy(sort(coeff_ncim(:,id)),(1:L3)/L3,'color',dgreen);
hold off
set(gca,'fontsize', fontsize);
xlabel('Coefficient value');
ylabel('CDF');

figure(23)
id1 = 1;
h1 = plot(coeff_mod_gravity(:,id1),'color',dgrey);
hold on
h2 = plot(coeff_ncim(:,id1),'color',dgreen,'linewidth',thick_line);
hold off
set(gca,'fontsize', fontsize);
xlim([0,n_TMs]);
xlabel('t_k');
ylabel('Coefficient magnitude');

if (print_figures)
    set(gcf, 'PaperPosition', paper_position);
    print(device, [outfile 'coeff_' num2str(id1) '_sigmasq_' num2str(sigmasq,'%.0e') '.eps']);
    fprintf('printing %s\n', [outfile 'coeff_' num2str(id1) '_sigmasq_' num2str(sigmasq,'%.0e') '.eps']);
end

figure(33)
id2 = 2;
h1 = plot(coeff_mod_gravity(:,id2),'color',dgrey);
hold on
h2 = plot(coeff_ncim(:,id2),'color',dgreen,'linewidth',thick_line);
hold off
set(gca,'fontsize', fontsize);
xlim([0,n_TMs]);
xlabel('t_k');
ylabel('Coefficient magnitude');

if (print_figures)
    set(gcf, 'PaperPosition', paper_position);
    print(device, [outfile 'coeff_' num2str(id2) '_sigmasq_' num2str(sigmasq,'%.0e') '.eps']);
    fprintf('printing %s\n', [outfile 'coeff_' num2str(id2) '_sigmasq_' num2str(sigmasq,'%.0e') '.eps']);
end

figure(43)
id3 = 143;
h1 = plot(coeff_mod_gravity(:,id3),'color',dgrey);
hold on
h2 = plot(coeff_ncim(:,id3),'color',dgreen,'linewidth',thick_line);
hold off
set(gca,'fontsize', fontsize);
xlim([0,n_TMs]);
xlabel('t_k');
ylabel('Coefficient magnitude');

if (print_figures)
    set(gcf, 'PaperPosition', paper_position);
    print(device, [outfile 'coeff_' num2str(id3) '_sigmasq_' num2str(sigmasq,'%.0e') '.eps']);
    fprintf('printing %s\n', [outfile 'coeff_' num2str(id3) '_sigmasq_' num2str(sigmasq,'%.0e') '.eps']);
end

%% Perform normal distribution fit on modulated gravity's coeff
figure(14)
id4 = 143;
L = length(coeff_mod_gravity(:,id4));
[muhat,sigmahat,muci,sigmaci] = normfit(coeff_mod_gravity(:,id4));
coeff_gaussian = normcdf(coeff_mod_gravity(:,id4),muhat,sigmahat);
[mod_grav_sorted,id_sort] = sort(coeff_mod_gravity(:,id4));

h1 = semilogy(mod_grav_sorted,coeff_gaussian(id_sort),'r','linewidth',thick_line);
hold on
h2 = semilogy(mod_grav_sorted,(1:L)/L,'color',dgrey,'linewidth',thick_line);
hold off
set(gca,'fontsize', fontsize);
xlabel('Coefficient value');
ylabel('CDF');
xlim([mod_grav_sorted(1) mod_grav_sorted(end)])
ylabel('CDF');

hlegend = legend([h1], 'Normal fit', 'Location', 'SouthEast');

if (print_figures)
    set(gcf, 'PaperPosition', paper_position);
    print(device, [outfile 'coeff_' num2str(id3) '_sigmasq_' num2str(sigmasq,'%.0e') '_normfit.eps']);
    fprintf('printing %s\n', [outfile 'coeff_' num2str(id3) '_sigmasq_' num2str(sigmasq,'%.0e') '_normfit.eps']);
end

%% Plot spatial coefficients
figure(15)
h1 = semilogy(sort(sqrt(spatial_mod_gravity),'descend'),'color',dgrey,'linewidth',thick_line);
hold on
h2 = semilogy(sort(sqrt(spatial_ncim),'descend'),'color',dgreen,'linewidth',thick_line);
hold off
set(gca,'fontsize', fontsize);
xlabel('Spatial coefficient order');
ylabel('Magnitude');
xlim([0 25]);
ylim([10^(-2) 10^(4)])

if (sigmasq == 1e-5)
    hlegend = legend([h1 h2], 'Modulated gravity', 'NCIM', 'Location','NorthEast');
end

if (print_figures)
    set(gcf, 'PaperPosition', paper_position);
    print(device, [outfile 'sigmasq_' num2str(sigmasq,'%.0e') '_spatial.eps']);
    fprintf('printing %s\n', [outfile 'sigmasq_' num2str(sigmasq,'%.0e') '_spatial.eps']);
end

figure(16)
h1 = semilogy(sort(sqrt(temporal_mod_gravity),'descend'),'color',dgrey,'linewidth',thick_line);
hold on
h2 = semilogy(sort(sqrt(temporal_ncim),'descend'),'color',dgreen,'linewidth',thick_line);
hold off
set(gca,'fontsize', fontsize);
xlabel('Temporal coefficient order');
ylabel('Magnitude');

max_x = max(sum(temporal_mod_gravity > 0),sum(temporal_ncim > 0));
xlim([0 25]);
ylim([10^(-2) 10^4]);

if (sigmasq == 1e-5)
    hlegend = legend([h1 h2], 'Modulated gravity', 'NCIM', 'Location','NorthEast');
end

if (print_figures)
    set(gcf, 'PaperPosition', paper_position);
    print(device, [outfile 'sigmasq_' num2str(sigmasq,'%.0e') '_temporal.eps']);
    fprintf('printing %s\n', [outfile 'sigmasq_' num2str(sigmasq,'%.0e') '_temporal.eps']);
end