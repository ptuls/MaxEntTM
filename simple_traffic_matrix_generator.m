%%SIMPLE_TRAFFIC_MATRIX_GENERATOR Generate a series of traffic matrices via
%   the Modulated Gravity Model.

% assume measurement interval is 15 minutes, so there are 24*4 = 96
% measurement intervals in one day

clear;
num_pops = 20;      % number of PoPs
day = 96;
num_days = 7;      % set number of days
n_TMs = num_days*day;        % number of traffic matrices
mean_traffic = 100; % average total traffic
pm_ratio = 2;       % peak-to-mean ratio
t_ratio = 0.25*pm_ratio;      % trough-to-mean ratio
diurnal_freq = 1/day;
spatial_var = 100;   % \sigma^2 parameter; variation of traffic
temporal_var = 0.01;
print_figures = 0;    % set to 1 to print figures

% generate total traffic
total_traffic = peak_mean_cycle(diurnal_freq,n_TMs,mean_traffic,pm_ratio,t_ratio);
if min(total_traffic) < 0
    % rescale
    total_traffic = total_traffic + abs(min(total_traffic));
    mean_traffic = mean(total_traffic);
end

% randomly generate incoming and outgoing total PoP traffic
% here we take the ratio of uniform random variables which turns out to be
% equivalent to a Beta distribution
% outgoing 
fraction = rand(num_pops,1);
fraction = fraction/(sum(fraction));

mean_row = fraction*mean_traffic;

% incoming
% fraction = rand(num_pops,1);
fraction = fraction/(sum(fraction));

mean_col = fraction*mean_traffic;

% note: rank of G must be 1
tic
[TM,G] = modulated_gravity(mean_row,mean_col,total_traffic,...
            spatial_var,temporal_var);
toc

% average gravity model: for user comparisons if desired
Gmean = mean_row*mean_col'/sum(mean_row);

%% Plot parameters
fontsize = 20;
linewidth = 1.1;
axeswidth = 3;
thickline = 2;
dgreen = [0.05 0.5 0.05];
dgrey = [0.65 0.65 0.65];

% total traffic from the generated traffic matrices
TM_total_traffic = squeeze(sum(sum(TM)));
% Note: the average total traffic over an interval will not be exactly
% equal to the constraint, since these are soft constraints. The average
% total traffic of the TMs would be determined by the total traffic of G
fprintf('\nStatistics\n');
fprintf('Average total traffic: %.2f\n',mean(TM_total_traffic));
fprintf('Standard deviation total traffic: %.2f\n',std(TM_total_traffic));
fprintf('Rank of gravity model: %d\n',rank(G));
fprintf('Largest principal component coefficient: %.2f\n', sum(diag(G)));

figure(1)
h1 = plot(total_traffic(1:day), 'color', dgrey, 'LineWidth', thickline);
hold on
h2 = plot(TM_total_traffic(1:day), 'color', 'r', 'LineWidth', linewidth);
hold off
xlim([0 day]);

hlegend = legend([h1 h2], 'Constraint', 'MGM', 'Location', 'NorthEast');
hx = xlabel('Interval');
hy = ylabel('Gbps');
axes_handle = gca;
set(axes_handle, 'FontSize', fontsize,'LineWidth',axeswidth);
set(hx, 'FontSize', fontsize);
set(hy, 'FontSize', fontsize);
% title('Total traffic in one day');

if print_figures
    device = '-depsc';
    suffix = 'eps';
    print(device, sprintf('mgm_sigmasq_%s_pm_%s.%s', ...
        num2str(spatial_var,'%d'),...
        num2str(pm_ratio,'%d'), suffix));
end