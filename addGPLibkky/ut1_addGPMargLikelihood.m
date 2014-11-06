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
f2 = @(x) -( 0.5*x.^3 + (x.^2) .*(sin(x)) + 1.5*(x-0.5).^2 );
f = @(x) f1(x(:,1)) + f2(x(:,2));
trueFuncHs = {f1, f2};
numGroups = 2;
bounds = [0 1; 0 1];

% Generate Data
m = 10;
X = rand(m, 2);
Y = f(X);

% Perform Additive GP Regression
decomposition = { [1], [2] };
hyperParams.useSameSm = false;
hyperParams.useSamePr = true;
hyperParams.fixSm = false;
hyperParams.fixPr = false;
hyperParams.sigmaSmRanges = [];
hyperParams.sigmaPrRange = [];
hyperParams.noises = 0.00 * std(Y) * ones(numGroups, 1);
hyperParams.commonNoise = 0.01 * std(Y);
hyperParams.meanFuncs = [];
hyperParams.commonMeanFunc = [];
dummyPt = zeros(0, 2);
[mu, KPost, Mus, KPosts, combinedFuncH, funcHs, sigmaSmOpts, sigmaPrOpts] = ...
  addGPMargLikelihood(X, Y, dummyPt, decomposition, hyperParams);
fprintf('Chosen Smoothness: %s\n', mat2str(sigmaSmOpts));
fprintf('Chosen Prior: %s\n', mat2str(sigmaPrOpts));
% Now plot the functions
plotAddGPOutputs(combinedFuncH, funcHs, f, trueFuncHs, decomposition, bounds);

% Now do some sanity checks on this - do a generic GP Regression
hp.meanFunc = [];
hp.sigmaSmRange = [];
hp.sigmaPrRange = [];
hp.sigmaSm = 0;
hp.sigmaPr = 0;
hp.noise = hyperParams.commonNoise;
hp.retFunc = true;
[~, ~, naiveGPFuncH] = GPMargLikelihood(X, Y, dummyPt, hp);

% Now do the test
numTestPts = 1000;
testPts = rand(numTestPts, 2);
trueVals = f(testPts);
addGPVals = combinedFuncH(testPts);
addGPErr = norm(trueVals-addGPVals)/sqrt(numTestPts); 
naiveGPVals = naiveGPFuncH(testPts);
naiveGPErr = norm(trueVals-naiveGPVals)/sqrt(numTestPts); 
fprintf('Add-GP Error: %0.5f\nNaive-GP Error:%0.5f\n', ...
  addGPErr, naiveGPErr);

