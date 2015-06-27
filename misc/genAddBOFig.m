% Generates some figures for BO/GPB
close all;
clear all;
addpath ~/libs/kky-matlab/utils/
addpath ~/libs/kky-matlab/ancillary/
addpath ../addGPLibkky/
rng(97);

N = 100;
numSamples = 20;
plotIdxs = [2:7 13:15 18:20];

func1 = @(t) - 70*(t-0).* (t-0.35).* (t+0.55).* (t-0.65).* (t-0.98);
func2 = @(t) - 70*(1.1-t).* (-t+0.6).* (-t+1.55).* (-t+0.4).* (-t+0.02);
func = @(t) func1(t(:,1)) + func2(t(:,2));

X = rand(20, 2);
Y = func(X);

c1 = [150 75 0]/255;
c2 = 'r';

LW = 3;
MS = 20;
LWS = 1.5;
FS = 22;

plotlims1 = [0 1 -0.6 1.50];
plotlims2 = [0 1 -0.6 1.90];
yLabPosn1 = [-0.08 1.40];
yLabPosn2 = [-0.08 1.80];
xLabPosn = [0.96 -0.70];
% figSize = [695 368 905 420];
% figSizeSq = [695 00 905 835];
figSize =   [751 380 835 420];
figSizeSq = [751 1 835 793];
saveFormat = 'eps';


% First plot the function
th = linspace(0,1,100)';
T = meshgrid(th, th);
figure; 
plot(th, func1(th), 'LineWidth', LW, 'Color', 'k'); hold on,
set(gca, 'Ytick', []);
axis(plotlims1);
xlabel('$x_{\{1\}}$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
% ylabel('$f^{(1)}(x_{\{1\}})$', 'rot', 0, 'Position', yLabPosn1, 'FontSize', FS, ...
%   'Interpreter', 'Latex');
text(0.2, 0.7, '$f^{(1)}(x_{\{1\}})$', 'FontSize', FS, ...
  'Interpreter', 'Latex');
set(gcf, 'Position', figSize);
% saveas(gcf, 'add1', saveFormat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(th, func2(th), 'LineWidth', LW, 'Color', 'k'); hold on,
set(gca, 'Ytick', []);
axis(plotlims2);
xlabel('$x_{\{2\}}$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
% ylabel('$f^{(2)}(x_{\{2\}})$', 'rot', 0, 'Position', yLabPosn2, 'FontSize', FS, ...
%   'Interpreter', 'Latex');
text(0.6, 1.0, '$f^{(2)}(x_{\{2\}})$', 'FontSize', FS, ...
  'Interpreter', 'Latex');
set(gcf, 'Position', figSize);
% saveas(gcf, 'add2', saveFormat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot2DFunction(func, [0 1; 0 1], 'contour'); hold on,
plot(X(plotIdxs, 1), X(plotIdxs, 2), 'kx', 'MarkerSize', MS, 'LineWidth', LW);
set(gcf, 'Position', figSizeSq);
set(gca, 'Ytick', []);
set(gca, 'Xtick', []);
set(gca, 'position', [0.01 0.01 0.98 0.98],'units','normalized')
% saveas(gcf, 'func2D', saveFormat);


decomposition = {[1], [2]};
hyperParams.useSameSm = true;
hyperParams.useSamePr = true;
hyperParams.fixSm = false;
hyperParams.fixPr = false;
hyperParams.sigmaSmRanges = [];
hyperParams.sigmaSmRange = [];
hyperParams.sigmaPrRange = [];
hyperParams.sigmaPrRanges = [];
hyperParams.noises = 0.00 * std(Y) * ones(2, 1);
hyperParams.commonNoise = 0.01 * std(Y);
hyperParams.meanFuncs = [];
hyperParams.commonMeanFunc = @(arg) mean(Y) * ones(size(arg,1), 1);
dummyPt = zeros(0, 2);
[mu, KPost, Mus, KPosts, combinedFuncH, funcHs, sigmaSmOpts, sigmaPrOpts] = ...
  addGPMargLikelihood(X, Y, dummyPt, decomposition, hyperParams);

[mu1, std1, K1] = funcHs{1}(th);
[mu2, std2, K2] = funcHs{2}(th);


figure;
gpSamples1 = GPDrawSamples(mu1, K1, numSamples);
gpSamples1 = gpSamples1 + 0.4;
gpSamples2 = GPDrawSamples(mu2, K2, numSamples);
gpSamples2 = gpSamples2 + 0.6;
% plot(th, func1(th), 'LineWidth', LW, 'Color', 'k'); hold on,
plot(th, gpSamples1, 'LineWidth', LWS);
axis(plotlims1);
xlabel('$x_{\{1\}}$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
set(gca, 'Ytick', []);
set(gcf, 'Position', figSize);
% saveas(gcf, 'gpPost1', saveFormat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
% plot(th, func2(th), 'LineWidth', LW, 'Color', 'k'); hold on,
plot(th, gpSamples2, 'LineWidth', LWS); hold on,
axis(plotlims2);
xlabel('$x_{\{2\}}$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
set(gcf, 'Position', figSize);
set(gca, 'Ytick', []);
% saveas(gcf, 'gpPost2', saveFormat);


% Plot UCB
figure;
betath = 2;

ucb1 = mu1 + betath * std1 + 0.3;
[maxVal1, maxIdx1] = max(ucb1); maxTh1 = th(maxIdx1);
plot(th, ucb1, 'Color', 'b', 'LineWidth', LW); hold on;
plot(maxTh1, maxVal1, 'r*', 'MarkerSize', MS, 'LineWidth', LW);
plot([maxTh1, maxTh1], [-0.7 maxVal1], 'r');
text(th(maxIdx1) - 0.25, maxVal1 - 0.01, '${\bf x^{(1)}_{t}} = 0.869$', ...
  'Color', 'r', 'FontSize', FS, 'Interpreter', 'Latex');
axis(plotlims1);
xlabel('$x_{\{1\}}$', 'Position', xLabPosn, 'FontSize', FS, ...
  'Interpreter', 'Latex');
% ylabel('$\tilde{\varphi}^{(1)}(x_{\{1\}})$', 'rot', 0, 'Position', ...
%   yLabPosn1, 'FontSize', FS, 'Interpreter', 'Latex');
text(0.20, 0.9, '$\tilde{\varphi}^{(1)}(x_{\{1\}})$', 'FontSize', FS, ...
  'Interpreter', 'Latex', 'Color', 'b');
set(gcf, 'Position', figSize);
set(gca, 'Ytick', []);
% saveas(gcf, 'ucb1', saveFormat);

figure;
ucb2 = mu2 + betath * std2 + 0.5;
[maxVal2, maxIdx2] = max(ucb2); maxTh2 = th(maxIdx2);
plot(th, ucb2, 'Color', 'b', 'LineWidth', LW); hold on;
plot(maxTh2, maxVal2, 'r*', 'MarkerSize', MS, 'LineWidth', LW);
plot([maxTh2, maxTh2], [-0.7 maxVal2], 'r');
text(th(maxIdx2) + 0.04, maxVal2 + 0.02, ...
  '${\bf x^{(2)}_{t}} = 0.141$', 'Interpreter', 'Latex', ...
  'FontSize', FS, 'Color', 'r');
axis(plotlims2);
xlabel('$x_{\{2\}}$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
% ylabel('$\tilde{\varphi}^{(2)}(x_{\{2\}})$', 'rot', 0, 'Position', ...
%   yLabPosn2, 'FontSize', FS, 'Interpreter', 'Latex');
text(0.6, 1.2, '$\tilde{\varphi}^{(2)}(x_{\{2\}})$', 'FontSize', FS, ...
  'Interpreter', 'Latex', 'Color', 'b');
set(gcf, 'Position', figSize);
set(gca, 'Ytick', []);
% saveas(gcf, 'ucb2', saveFormat);

figure;
plot2DFunction(func, [0 1; 0 1], 'contour'); hold on,
plot(X(plotIdxs, 1), X(plotIdxs, 2), 'kx', 'MarkerSize', MS, 'LineWidth', LW);
plot(maxTh1, maxTh2, 'r*', 'MarkerSize', MS, 'LineWidth', LW);
text(0.6, 0.175, '${\bf{x_{t}}} = (0.869,0.141)$', 'Interpreter', 'Latex', ...
  'FontSize', FS, 'Color', 'r');
set(gcf, 'Position', figSizeSq);
set(gca, 'Ytick', []);
set(gca, 'Xtick', []);
set(gca,'position',[0.01 0.01 0.98 0.98],'units','normalized')
% saveas(gcf, 'pts2Dnext', saveFormat);

figure;
plot2DFunction(func, [0 1; 0 1], 'contour'); hold on,
set(gcf, 'Position', figSizeSq);
set(gca, 'Ytick', []);
set(gca, 'Xtick', []);
set(gca,'position',[0.01 0.01 0.98 0.98],'units','normalized')
% saveas(gcf, 'func2D', saveFormat);
