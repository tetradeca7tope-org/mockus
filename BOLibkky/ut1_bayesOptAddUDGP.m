% Unit Test for additive GP Bayesian Optimization with unknown decomposition
addpath ~/libs/kky-matlab/GPLibkky/
addpath ~/libs/kky-matlab/utils/
addpath ~/libs/kky-matlab/ancillary/
addpath ../addGPLibkky/
addpath ../utils/

close all;
clear all;
clc; clc;

% Test 1
% 2 1D examples
fprintf('Test 1\n======================\n');
f1 = @(x) x.^2 + x.*sin(4*x);
f2 = @(x) -( 0.5*x.^3 + (x.^2) .*(sin(x)) + 1.5*(x-0.5).^2 );
func = @(x) f1(x(:,1)) + f2(x(:,2));
trueFuncHs = {f1, f2};
bounds = [0 1; 0 1];

% Prelims
numGroups = 2;
decomposition = { [1], [2] };
numDims = 2;
numIters = 40;

% Set the parameters for additive Bayesian Optimization
boAddParams.utilityFunc = 'UCB';
boAddParams.optPtStdThreshold = 0.002;
boAddParams.numInitPts = 4 * numDims;
boAddParams.meanFuncs = [];
boAddParams.commonMeanFunc = [];
boAddParams.commonNoise = [];
boAddParams.noises = 0 * ones(numGroups, 1);
boAddParams.fixPr = false;
boAddParams.useSamePr = true;
boAddParams.sigmaPrRange = []; % Let the Marginal Likelihood pick this
boAddParams.useFixedBandWidth = false;
boAddParams.fixSm = false;
boAddParams.useSameSm = true;
boAddParams.alBWLB = 1e-5;
boAddParams.alBWUB = 5;

% Call BO
[maxVal, maxPt, boAddQueries, boAddVals, boAddHist] = ...
  bayesOptAddUDGP(func, 1, 2, bounds, numIters, boAddParams);

% Plot the iteration statistics
figure;
plot( boAddHist, 'b');

% Plot the maxPt
figure;
plot2DFunction( func, bounds, 'contour'); hold on,
plot(maxPt(1), maxPt(2), 'rx', 'MarkerSize', 10, 'LineWidth', 3);
title('func');

