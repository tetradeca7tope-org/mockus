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
MS = 15;
LWS = 1.5;

c1 = [150 75 0]/255;
c2 = 'r';


% First plot the function
th = linspace(0,1,100)';
plot(th, func(th), 'LineWidth', LW, 'Color', 'k'); hold on,

% X = [0.1 0.3 0.6 0.77 0.95]';
% X = [0.1 0.4 0.75 0.95]';
X = [0.1 0.4 0.55 0.78 0.95]';
% X = [0.1 0.4 0.70 0.80 0.95]';
% X = [0.1 0.4 0.6 0.8 0.98]';
Y = func(X);
plot(X, Y, 'kx', 'MarkerSize', MS, 'LineWidth', LW);

% set(gca, 'Xtick', []);
set(gca, 'Ytick', []);

hyperparams.meanFunc = [];
hyperparams.sigmaSmRange = [];
hyperparams.sigmaPrRange = [];
hyperparams.sigmaSm = 0;
hyperparams.sigmaPr = 0;
hyperparams.noise = 0.0001;

figure;
plot(th, func(th), 'LineWidth', LW, 'Color', 'k'); hold on,
plot(X, Y, 'kx', 'MarkerSize', MS, 'LineWidth', LW);
[mu, K, funcH] = GPMargLikelihood(X, Y, th, hyperparams);
gpSamples = GPDrawSamples(mu, K, numSamples);
plot(th, gpSamples, 'LineWidth', LWS);

% set(gca, 'Xtick', []);
set(gca, 'Ytick', []);
% set(gca,'position',[0 0 1 1],'units','normalized');

% Plot UCB
figure;
betath = 2;
ucb = mu + betath * diag(K);
plot(th, ucb, 'Color', 'b', 'LineWidth', LW);
hold on;
[maxVal, maxIdx] = max(ucb);
plot(th(maxIdx), maxVal, 'r*', 'MarkerSize', MS, 'LineWidth', LW);
plot([th(maxIdx), th(maxIdx)], [-0.5 maxVal], 'r');
th(maxIdx),

text(th(maxIdx) - 0.23, maxVal + 0.025, ...
'\fontsize{20} \color{red} {\bf x_{t+1}} = 0.828');

% set(gca, 'Xtick', []);
set(gca, 'Ytick', []);
% set(gca,'position',[0 0 1 1],'units','normalized');
