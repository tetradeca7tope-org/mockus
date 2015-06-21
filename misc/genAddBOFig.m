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


LW = 3;
MS = 15;
LWS = 1.5;

c1 = [150 75 0]/255;
c2 = 'r';


% First plot the function
th = linspace(0,1,100)';
T = meshgrid(th, th);
figure; 
plot(th, func1(th), 'LineWidth', LW, 'Color', 'k'); hold on,
set(gca, 'Ytick', []);
figure;
plot(th, func2(th), 'LineWidth', LW, 'Color', 'k'); hold on,
set(gca, 'Ytick', []);
figure;
plot2DFunction(func, [0 1; 0 1], 'contour'); hold on,
plot(X(plotIdxs, 1), X(plotIdxs, 2), 'kx', 'MarkerSize', MS, 'LineWidth', LW);
set(gca, 'Ytick', []);
set(gca, 'Xtick', []);


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
gpSamples2 = GPDrawSamples(mu2, K2, numSamples);
plot(th, gpSamples1, 'LineWidth', LWS);
set(gca, 'Ytick', []);
figure;
plot(th, gpSamples2, 'LineWidth', LWS);
set(gca, 'Ytick', []);

% set(gca, 'Xtick', []);
set(gca, 'Ytick', []);
% set(gca,'position',[0 0 1 1],'units','normalized');

% Plot UCB
figure;
betath = 2;

ucb1 = mu1 + betath * std1;
[maxVal1, maxIdx1] = max(ucb1); maxTh1 = th(maxIdx1);
plot(th, ucb1, 'Color', 'b', 'LineWidth', LW); hold on;
plot(maxTh1, maxVal1, 'r*', 'MarkerSize', MS, 'LineWidth', LW);
plot([maxTh1, maxTh1], [-0.7 maxVal1], 'r');
text(th(maxIdx1) - 0.23, maxVal1 - 0.12, ...
'\fontsize{20} \color{red} {\bf x^{(1)}_{t+1}} = 0.869');
ylim([-0.6 1.2]);
set(gca, 'Ytick', []);

figure;
ucb2 = mu2 + betath * std2;
[maxVal2, maxIdx2] = max(ucb2); maxTh2 = th(maxIdx2);
plot(th, ucb2, 'Color', 'b', 'LineWidth', LW); hold on;
plot(maxTh2, maxVal2, 'r*', 'MarkerSize', MS, 'LineWidth', LW);
plot([maxTh2, maxTh2], [-0.7 maxVal2], 'r');
text(th(maxIdx2) + 0.03, maxVal2 - 0.12, ...
'\fontsize{20} \color{red} {\bf x^{(2)}_{t+1}} = 0.141');
ylim([-0.6 1.3]);
set(gca, 'Ytick', []);

figure;
plot2DFunction(func, [0 1; 0 1], 'contour'); hold on,
plot(X(plotIdxs, 1), X(plotIdxs, 2), 'kx', 'MarkerSize', MS, 'LineWidth', LW);
plot(maxTh1, maxTh2, 'r*', 'MarkerSize', MS, 'LineWidth', LW);
text(0.4, 0.18, ...
'\fontsize{20} \color{red} {\bf x_{t+1}} = (0.869,0.141)');
set(gca, 'Ytick', []);
set(gca, 'Xtick', []);

figure;
plot2DFunction(func, [0 1; 0 1], 'contour'); hold on,
set(gca, 'Ytick', []);
set(gca, 'Xtick', []);
