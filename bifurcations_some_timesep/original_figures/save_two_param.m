var1 = 'varphi';
var2 = 'beta_G';

set(0,'defaultaxesfontsize', 25);
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesLineWidth',2);
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on');
set(0,'defaultLineMarkerSize', 25);
set(0, 'DefaultFigureRenderer', 'painters')


xlabel(['$\' var1 '$'], 'Interpreter', 'latex')
ylabel(['$\' var2 '$'], 'Interpreter', 'latex')

saveas(gcf, [var1 '_' var2 '.png']);
saveas(gcf, [var1 '_' var2 '.svg']);