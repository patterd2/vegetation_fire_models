% Plot proportions of FGBA (in data variable)
variable = 'varphi';
folder = 'bifurcations_some_timesep';

% Default settings
set(0, 'DefaultFigureRenderer', 'painters');
set(0, 'defaultfigureposition', [400 250 500 500])
set(0,'defaultaxesfontsize', 25);
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesLineWidth',2);
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')
set(0,'defaultlineMarkerSize', 25)

% Figure
close(figure(1))
figure(1)

% Append bottom corners to data
var_min = min(data(1,:));
var_max = max(data(1,:));
data_new = horzcat([var_min, 0, 0, 1]', data, [var_max, 0, 0, 1]');

% Compute FGBA lines
num_points = size(data_new, 2);
F = data_new(2, :);
GF = data_new(2, :) + data_new(3, :);
AGF = 1 - data_new(4, :);

% Set up color map
custom_map = [0.4660 0.6740 0.1880
    0.0039 0.1953 0.1250
    0.8500 0.3250 0.0980
    0.5 0.5 0.5]; %light green, dark green, orange, gray colors

fill([var_min, var_max, var_max, var_min], [0, 0, 1, 1], [0.8500 0.3250 0.0980])
hold on
fill(data_new(1, :), AGF, [0.5 0.5 0.5])
fill(data_new(1, :), GF, [0.4660 0.6740 0.1880])
fill(data_new(1, :), F, [0.0039 0.1953 0.1250])

set(figure(1), 'position', [400 250 500 400])
xlabel(['$\' variable '$'], 'Interpreter', 'latex')
ylabel('Land fraction')
saveas(gcf, [folder '/' variable '_FGBA.png']);
saveas(gcf, [folder '/' variable '_FGBA.svg']);
save([folder '/' variable '.mat'], 'data','-v7.3')
