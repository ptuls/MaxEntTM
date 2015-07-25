%%ABILENE_TRAFFIC_MATRIX_EXAMPLE Generate a series of traffic matrices via
%   the Modulated Gravity Model using total traffic obtained by extracting
%   a normalized Minimum Description Length (nMDL) selection of frequency
%   coefficients from Abilene traffic, starting 

% a measurement interval is 5 minutes, so there are 24*12 = 288
% measurement intervals in one day

clear;
num_pops = 20;      % number of PoPs
spatial_var = 10;   % \sigma^2 parameter; variation of traffic
temporal_var = 1;
print_figures = 0;    % set to 1 to print figures

% generate total traffic
load data/abilene_data
total_traffic = abilene_traffic_nmdl;
mean_traffic = mean(total_traffic);
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
fontsize = 22;
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

day = 24*60/(5);

figure(1)
h1 = plot(total_traffic(1:7*day), 'color', dgrey, 'LineWidth', thickline);
hold on
h2 = plot(TM_total_traffic(1:7*day), 'color', 'r', 'LineWidth', linewidth);
hold off

% hlegend = legend([h1 h2], 'Constraint', 'MGM', 'Location', 'NorthWest');
hx = xlabel('t_k');
hy = ylabel('Gbps');
axes_handle = gca;
set(axes_handle, 'FontSize', fontsize,'LineWidth',axeswidth);
set(hx, 'FontSize', fontsize);
set(hy, 'FontSize', fontsize);
set(gca,'fontsize', fontsize, 'XTick', day*(0:6), ...
    'XTickLabel', {'Mon','Tue','Wed','Thu','Fri','Sat','Sun'});
xlim([0 2000]);
xlabel('Week starting 2004-03-01');
ylabel('GB');

if print_figures
    paper_position = [0 0 14 12];
    set(gcf, 'PaperPosition', paper_position);
    device = '-depsc';
    suffix = 'eps';
    print(device, sprintf('abilene_mgm_sigmasq_%s_pm_%s.%s', ...
        num2str(spatial_var,'%d'),...
        num2str(pm_ratio,'%f'), suffix));
end