% Unit Test for additive GP Regression

addpath ~/libs/kky-matlab/GPLibkky/
addpath ~/libs/kky-matlab/utils/
addpath ~/libs/kky-matlab/ancillary/
addpath ../BOLibkky/
addpath ../utils/

close all;
clear all;

% Test 1
% 2 1D examples
fprintf('Test 1\n======================\n');
f1 = @(x) x.^2 + x.*sin(4*x);
f2 = @(x) x.^3 + (x.^2) .* sin(x);
f = @(x) f1(x(:,1)) + f2(x(:,2));
trueFuncHs = {f1, f2};
bounds = [0 1; 0 1];

% Generate Data
m = 20;
X = rand(m, 2);
Y = f(X);

% Perform Additive GP Regression
decomposition = { [1], [2] };
hyperParams.sigmaSms = 0.3 * norm( std(X) );
hyperParams.sigmaPrs = 2 * std(Y);
hyperParams.noises = 0.05 * std(Y);
hyperParams.commonNoise = 0.1 * std(Y);
hyperParams.meanFuncs = [];
hyperParams.commonMeanFunc = [];
dummyPt = zeros(0, 2);
[mu, KPost, Mus, KPosts, combinedFuncH, funcHs] = ...
  addGPRegression(X, Y, dummyPt, decomposition, hyperParams);
% Now plot the functions
plotAddGPOutputs(combinedFuncH, funcHs, f, trueFuncHs, decomposition, bounds);

