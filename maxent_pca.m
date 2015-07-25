%% MAXENT_PCA Compute the principle component analysis of the
%  Modulated Gravity Model (MGM) and the Non-stationary Conditionally
%  Independent Model (NCIM). 

%
% file:     maxent_pca.m, (c) Paul Tune, Mon Dec 22 2014
% created: 	Mon Dec 22 2014 
% author:  	Paul Tune 
% email:   	paul.tune@adelaide.edu.au

clear;

% load Abilene traffic matrices; we will extract the Fourier coefficients
% from these series of traffic matrices
load tm_real

% turn into 3D array and compute row and column marginals
tm_tensor = reshape(tm_real,N,N,n_TMs);
r = squeeze(sum(tm_tensor,2));  % row marginals
c = squeeze(sum(tm_tensor,1));  % column marginals

% total traffic in each time interval
T = sum(tm_real,1);

% compute mean
mean_r = mean(r,2);
mean_c = mean(c,2);
mean_T = mean(T);

% extract Fourier coefficients and form k-sparse version
fT = fft(T-mean_T);
[fTmag,indf] = sort(abs(fT(1:(n_TMs/2))),2,'descend');

% choose k coefficients based on nMDL
k = 23;
fs = zeros(1,n_TMs/2);
fs(indf(1:k)) = fT(indf(1:k));
fsparse = [fs fT(n_TMs/2+1) conj(fs(end:-1:2))];
norm_T = real(ifft(fsparse));

%% noise parameter \sigma^2 for both models
sigmasq = 10;

%% modulated gravity based on sparse Fourier, nMDL selected coefficients
rate = 5;

X = zeros(num_nodes,num_nodes,n_TMs);
Tk = mean(T);
rk = mean(r,2);
ck = mean(c,2);




tic
pU = rk/Tk;
pV = ck/Tk;
U = HMC_exact(eye(N), zeros(N,1), eye(N)*sigmasq/Tk^2, pU, true, 2, pU);
V = HMC_exact(eye(N), zeros(N,1), eye(N)*sigmasq/Tk^2, pV, true, 2, pV);
M0 = (Tk+norm_T)/Tk;
M = HMC_exact(eye(n_TMs),zeros(n_TMs,1),sigmasq/Tk^2,M0',true,2,M0');

G = Tk*U(:,2)*V(:,2)';
X = reshape(reshape(G,N^2,1)*M(:,2)',N,N,n_TMs);
toc

%% non-stationary conditionally independent model (NCIM) with delayed rows,
%  columns to represent time zone shifts

% need to construct the row and column spatiotemporal signals
% choose New York's as reference
idNY = 9;
rs = zeros(N,n_TMs);
cs = zeros(N,n_TMs);
rs(idNY,:) = freq_sparsify(r(idNY,:),k); % nMDL calculated optimal choice of frequencies
cs(idNY,:) = freq_sparsify(c(idNY,:),k);

% each measurement interval is 5 minutes, so need to convert hours to
% minutes
m = 60/5;

% delay due to time zones (everywhere except NY)
D = m*[0 0 1 2 1 0 1 3 3 3 0];
location = [1:8 10:N]; % everywhere except NY

% build constraints
rscale = mean_r/mean_r(idNY);
cscale = mean_c/mean_c(idNY);
for i=1:length(location)
    rs(location(i),:) = circshift(rscale(location(i))*rs(idNY,:),[0 D(i)]);
    cs(location(i),:) = circshift(cscale(location(i))*cs(idNY,:),[0 D(i)]);
end
Ts = sum(rs,1);

tic
for i=1:n_TMs 
    pU = rs(:,i)/Ts(i);
    pV = cs(:,i)/Ts(i);
    
    U = HMC_exact(eye(N), zeros(N,1), eye(N)*sigmasq, pU, true, 2, pU);
    V = HMC_exact(eye(N), zeros(N,1), eye(N)*sigmasq, pV, true, 2, pV);
    Y(:,i) = U(:,2);
    Z(:,i) = V(:,2);
    Yd(:,i) = pU;
    Zd(:,i) = pV;
    Xc(:,:,i) =  Ts(i)*U(:,2)*V(:,2)';
end
toc

%% Compute traffic decomposition via PCA
x_real = reshape(tm_real,N^2,n_TMs);
x_mod_gravity = reshape(X,N^2,n_TMs);
x_nscim = reshape(Xc,N^2,n_TMs);
% x_g = reshape(Xg,N^2,n_TMs);

mean_x_real = mean(x_real,2);
mean_x_mod_gravity = mean(x_mod_gravity,2);
mean_x_nscim = mean(x_nscim,2);
% mean_x_g = mean(x_g,2);

coeff_real = princomp(x_real);
coeff_mod_gravity = princomp(x_mod_gravity);
coeff_nscim = princomp(x_nscim);

spatial_real = svd(x_real);
[spatial_mod_gravity] = eig(x_mod_gravity*x_mod_gravity'-mean_x_mod_gravity *mean_x_mod_gravity');
spatial_nscim = eig(x_nscim*x_nscim'-mean_x_nscim*mean_x_nscim');
% spatial_g = eig(x_g*x_g'-mean_x_g*mean_x_g');

temporal_real = eig(x_real'*x_real);
temporal_mod_gravity = eig(x_mod_gravity'*x_mod_gravity - mean_x_mod_gravity'*mean_x_mod_gravity);
temporal_nscim = eig(x_nscim'*x_nscim - mean_x_nscim'*mean_x_nscim);
% temporal_g = eig(x_g'*x_g - mean_x_g'*mean_x_g);


% remove machine noise
spatial_real(spatial_real < 0) = 0;
spatial_mod_gravity(spatial_mod_gravity < 0) = 0;
spatial_nscim(spatial_nscim < 0) = 0;

% remove machine noise
temporal_real(temporal_real < 0) = 0;
temporal_mod_gravity(temporal_mod_gravity < 0) = 0;
temporal_nscim(temporal_nscim < 0) = 0;

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
plot(squeeze(sum(sum(X))),'color',dgrey);
plot(squeeze(sum(sum(Xc))),'color',dgreen);
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
plot(coeff_nscim(:,1),'color',dgreen);
hold off
set(gca,'fontsize', fontsize);
xlim([0,n_TMs]);
xlabel('t_k');
ylabel('Coefficient magnitude');

figure(11)
plot(coeff_real(:,2),'color',dgrey,'linewidth',thick_line);
hold on
plot(coeff_mod_gravity(:,2),'color',dgrey);
plot(coeff_nscim(:,2),'color',dgreen);
hold off
set(gca,'fontsize', fontsize);
xlim([0,n_TMs]);
xlabel('t_k');
ylabel('Coefficient magnitude');

figure(12)
plot(coeff_real(:,3),'color',dgrey,'linewidth',thick_line);
hold on
plot(coeff_mod_gravity(:,3),'color',dgrey);
plot(coeff_nscim(:,3),'color',dgreen);
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
L3 = length(coeff_nscim(:,id));

semilogy(sort(coeff_real(:,id)),(1:L1)/L1,'color',dgrey,'linewidth',thick_line);
hold on
semilogy(sort(coeff_mod_gravity(:,id)),(1:L2)/L2,'r');
semilogy(sort(coeff_nscim(:,id)),(1:L3)/L3,'color',dgreen);
hold off
set(gca,'fontsize', fontsize);
xlabel('Coefficient value');
ylabel('CDF');

figure(23)
id1 = 1;
h1 = plot(coeff_mod_gravity(:,id1),'color',dgrey);
hold on
h2 = plot(coeff_nscim(:,id1),'color',dgreen,'linewidth',thick_line);
hold off
set(gca,'fontsize', fontsize);
xlim([0,n_TMs]);
xlabel('t_k');
ylabel('Coefficient magnitude');
ylim([0.01 0.035]);
if (sigmasq == 1e-5)
    hlegend = legend([h1 h2], 'Modulated gravity', 'NCIM', 'Location','NorthWest');
end

if (print_figures)
    set(gcf, 'PaperPosition', paper_position);
    print(device, [outfile 'coeff_' num2str(id1) '_sigmasq_' num2str(sigmasq,'%.0e') '.eps']);
    fprintf('printing %s\n', [outfile 'coeff_' num2str(id1) '_sigmasq_' num2str(sigmasq,'%.0e') '.eps']);
end

figure(33)
id2 = 2;
h1 = plot(coeff_mod_gravity(:,id2),'color',dgrey);
hold on
h2 = plot(coeff_nscim(:,id2),'color',dgreen,'linewidth',thick_line);
hold off
set(gca,'fontsize', fontsize);
xlim([0,n_TMs]);
xlabel('t_k');
ylabel('Coefficient magnitude');
ylim([0.01 0.035]);
% hlegend = legend([h1 h2], 'Modulated gravity', 'NCIM');

if (print_figures)
    set(gcf, 'PaperPosition', paper_position);
    print(device, [outfile 'coeff_' num2str(id2) '_sigmasq_' num2str(sigmasq,'%.0e') '.eps']);
    fprintf('printing %s\n', [outfile 'coeff_' num2str(id2) '_sigmasq_' num2str(sigmasq,'%.0e') '.eps']);
end

figure(43)
id3 = 143;
h1 = plot(coeff_mod_gravity(:,id3),'color',dgrey);
hold on
h2 = plot(coeff_nscim(:,id3),'color',dgreen,'linewidth',thick_line);
hold off
set(gca,'fontsize', fontsize);
xlim([0,n_TMs]);
xlabel('t_k');
ylabel('Coefficient magnitude');
ylim([-0.1 0.1]);
% hlegend = legend([h1 h2], 'Modulated gravity', 'NCIM');

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
% ylim([-0.1 0.1]);

if (print_figures)
    set(gcf, 'PaperPosition', paper_position);
    print(device, [outfile 'coeff_' num2str(id3) '_sigmasq_' num2str(sigmasq,'%.0e') '_normfit.eps']);
    fprintf('printing %s\n', [outfile 'coeff_' num2str(id3) '_sigmasq_' num2str(sigmasq,'%.0e') '_normfit.eps']);
end

figure(15)
h1 = semilogy(sort(sqrt(spatial_mod_gravity),'descend'),'color',dgrey,'linewidth',thick_line);
hold on
h2 = semilogy(sort(sqrt(spatial_nscim),'descend'),'color',dgreen,'linewidth',thick_line);
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
h2 = semilogy(sort(sqrt(temporal_nscim),'descend'),'color',dgreen,'linewidth',thick_line);
hold off
set(gca,'fontsize', fontsize);
xlabel('Temporal coefficient order');
ylabel('Magnitude');

max_x = max(sum(temporal_mod_gravity > 0),sum(temporal_nscim > 0));
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
% 
% figure(17)
% h1 = semilogy(sort(spatial_mod_gravity,'descend'),'color',dgrey,'linewidth',thick_line);
% hold on
% h2 = semilogy(sort(spatial_g,'descend'),'color','b','linewidth',thick_line);
% hold off
% set(gca,'fontsize', fontsize);
% xlabel('Spatial coefficient order');
% ylabel('Magnitude');
% 
% max_x = max(sum(temporal_mod_gravity > 0),sum(temporal_nscim > 0));
% % xlim([0 25]);
% % ylim([10^3 10^7]);
% 
% if (sigmasq == 1e-5)
%     hlegend = legend([h1 h2], 'Modulated gravity', 'Outer', 'Location','NorthEast');
% end
% 
% if (print_figures)
%     set(gcf, 'PaperPosition', paper_position);
%     print(device, [outfile 'sigmasq_' num2str(sigmasq,'%.0e') '_outer_spatial.eps']);
%     fprintf('printing %s\n', [outfile 'sigmasq_' num2str(sigmasq,'%.0e') '_outer_spatial.eps']);
% end
% 
% 
% figure(18)
% h1 = semilogy(sort(spatial_mod_gravity,'descend'),'color',dgrey,'linewidth',thick_line);
% hold on
% h2 = semilogy(sort(spatial_g,'descend'),'color','b','linewidth',thick_line);
% hold off
% set(gca,'fontsize', fontsize);
% xlabel('Temporal coefficient order');
% ylabel('Magnitude');
% 
% max_x = max(sum(temporal_mod_gravity > 0),sum(temporal_nscim > 0));
% % xlim([0 25]);
% % ylim([10^3 10^7]);
% 
% if (sigmasq == 1e-5)
%     hlegend = legend([h1 h2], 'Modulated gravity', 'Outer', 'Location','NorthEast');
% end
% 
% if (print_figures)
%     set(gcf, 'PaperPosition', paper_position);
%     print(device, [outfile 'sigmasq_' num2str(sigmasq,'%.0e') '_outer_temporal.eps']);
%     fprintf('printing %s\n', [outfile 'sigmasq_' num2str(sigmasq,'%.0e') '_outer_temporal.eps']);
% end
