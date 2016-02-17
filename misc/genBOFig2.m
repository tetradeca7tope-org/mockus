% Generates some figures for BO/GPB
close all;
clear all;
addpath ~/libs/kky-matlab/GPLibkky/
addpath ~/libs/kky-matlab/utils/
rng(104);

N = 100;
numSamples = 20;

func = @(t) - 70*(t-0).* (t-0.35).* (t+0.55).* (t-0.65).* (t-0.98);


LW = 3;
MS = 25;
LWS = 1.5;
FS = 26;

plotlims = [0 1 -0.6 1.30];
yLabPosn = [-0.05 1.20];
xLabPosn = [0.95 -0.68];
figSize = [85 261 898 420];
saveFormat = 'eps';

c1 = [150 75 0]/255;
c2 = 'r';

% First plot the function
th = linspace(0,1,100)';
[maxVal, maxIdx] = max(func(th));
maxPt = th(maxIdx);
plot(th, func(th), 'LineWidth', LW, 'Color', 'k'); hold on,
plot(maxPt, maxVal, 'r*', 'MarkerSize', MS, 'LineWidth', 1.5*LW);
plot([maxPt, maxPt], [-0.6 maxVal], 'r--', 'LineWidth', 2);
set(gca, 'Ytick', []);
set(gca, 'Xtick', []);
text(maxPt - 0.05, -0.52, '$x_*$', 'FontSize', FS, ...
  'Interpreter', 'Latex');
text(maxPt - 0.12, maxVal + 0.05, '$f(x_*)$', 'FontSize', 22,  ...
  'Interpreter', 'Latex');
axis(plotlims);
xlabel('$x$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
ylabel('$f(x)$', 'rot', 0, 'Position', yLabPosn, 'FontSize', FS, ...
  'Interpreter', 'Latex');
set(gcf, 'Position', figSize);
% saveas(gcf, 'func', saveFormat);
% ylabh = get(gca, 'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') - [0.2 0])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First plot the function
figure;
th = linspace(0,1,100)';
plot(th, func(th), 'LineWidth', LW, 'Color', 'k'); hold on,

% X = [0.1 0.3 0.6 0.77 0.95]';
% X = [0.1 0.4 0.75 0.95]';
X = [0.1 0.4 0.55 0.78 0.95]';
% X = [0.1 0.4 0.70 0.80 0.95]';
% X = [0.1 0.4 0.6 0.8 0.98]';
Y = func(X);
plot(X, Y, 'kx', 'MarkerSize', MS, 'LineWidth', 1.5*LW);
axis(plotlims);
xlabel('$x$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
ylabel('$f(x)$', 'rot', 0, 'Position', yLabPosn, 'FontSize', FS, ...
  'Interpreter', 'Latex');

% set(gca, 'Xtick', []);
set(gca, 'Ytick', []);
set(gcf, 'Position', figSize);
% saveas(gcf, 'funcPts', saveFormat);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hyperparams.meanFunc = [];
hyperparams.sigmaSmRange = [];
hyperparams.sigmaPrRange = [];
hyperparams.sigmaSm = 0;
hyperparams.sigmaPr = 0;
hyperparams.noise = 0.0001;

figure;
plot(th, func(th), 'LineWidth', LW, 'Color', 'k'); hold on,
[mu, K, funcH] = GPMargLikelihood(X, Y, th, hyperparams);
gpSamples = GPDrawSamples(mu, K, numSamples);
plot(th, gpSamples, 'LineWidth', LWS);
plot(X, Y, 'kx', 'MarkerSize', MS, 'LineWidth', 1.5*LW);

% set(gca, 'Xtick', []);
set(gca, 'Ytick', []);
axis(plotlims);
% set(gca,'position',[0 0 1 1],'units','normalized');
xlabel('$x$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
ylabel('$f(x)$', 'rot', 0, 'Position', yLabPosn, 'FontSize', FS, ...
  'Interpreter', 'Latex');
set(gcf, 'Position', figSize);
% saveas(gcf, 'gpPost', saveFormat);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot UCB
figure;
betath = 2;
ucb = mu + betath * diag(K);
plot(th, ucb, 'Color', 'b', 'LineWidth', LW);
hold on;
[maxVal, maxIdx] = max(ucb);
plot(th(maxIdx), maxVal, 'r*', 'MarkerSize', MS, 'LineWidth', 1.5*LW);
plot([th(maxIdx), th(maxIdx)], [-0.5 maxVal], 'r');
th(maxIdx),

text(0.6, -0.32, ...
'${\bf{x_{t}}} = 0.828$', 'Color', 'r', 'Interpreter', 'Latex', ...
  'FontSize', FS);

% set(gca, 'Xtick', []);
set(gca, 'Ytick', []);
axis([0 1 -0.5 1.05]);
% set(gca,'position',[0 0 1 1],'units','normalized');
xlabel('$x$', 'Position', [0.95 -0.60], 'FontSize', FS, 'Interpreter','Latex');
ylabel('$\varphi_t(x)$', 'Interpreter', 'Latex', ...
   'rot', 0, 'Position', [-0.05 1], 'FontSize', FS);
set(gcf, 'Position', figSize);
% saveas(gcf, 'ucb', saveFormat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Function and UCB together
figure;
% X = [0.1 0.3 0.6 0.77 0.95]';
% X = [0.05 0.4 0.70 0.95]';
X = [0.05 0.3 0.5 0.67 0.95]';
% X = [0.03 0.3 0.45 0.7 0.95]';
% X = [0.1 0.4 0.70 0.80 0.95]';
% X = [0.1 0.4 0.6 0.8 0.98]';
Y = func(X);
[mu, K, funcH] = GPMargLikelihood(X, Y, th, hyperparams);

FS = 50;

betath = 4;
ucb = mu + betath * sqrt(diag(K));
[maxVal, maxIdx] = max(ucb);
plot(th, func(th), 'LineWidth', LW, 'Color', 'k'); hold on,
plot(th, ucb, '--b', 'LineWidth', LW);
plot(th(maxIdx), maxVal, 'r*', 'MarkerSize', MS, 'LineWidth', 1.5*LW);

Y = func(X);
plot(X, Y, 'kx', 'MarkerSize', MS, 'LineWidth', 1.5*LW);
newPlotLims = [0 1 -0.6 1.7];
axis(newPlotLims);
xlabel('$x$', 'Position', [0.95 -0.4], 'FontSize', FS, 'Interpreter', 'Latex');
% ylabel('$f(x)$', 'rot', 0, 'Position', [-0.05 1.44], 'FontSize', FS, ...
%   'Interpreter', 'Latex');
text(0.215, 1.44, '$\varphi_t$', 'FontSize', FS, ...
  'Interpreter', 'Latex', 'Color', 'b');
text(0.16, 0.32, '$f$', 'FontSize', FS, ...
  'Interpreter', 'Latex');

% plot([th(maxIdx), th(maxIdx)], [-0.5 maxVal], 'r');
% text(0.8, -0.32, ...
% '${\bf{x_{t}}}$', 'Color', 'r', 'Interpreter', 'Latex', ...
%   'FontSize', FS);

set(gca, 'Xtick', []);
set(gca, 'Ytick', []);
set(gcf, 'Position', figSize);
% saveas(gcf, 'funcPts', saveFormat);
    set(gca, 'Position', [0.02 0.02 0.96 0.96], 'units', 'normalized');
    set(gcf, 'Position', figSize);



