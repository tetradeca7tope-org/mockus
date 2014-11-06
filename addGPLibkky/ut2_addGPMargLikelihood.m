% Unit Test for additive GP Regression

addpath ~/libs/kky-matlab/GPLibkky/
addpath ~/libs/kky-matlab/utils/
addpath ~/libs/kky-matlab/ancillary/
addpath ../BOLibkky/
addpath ../utils/

close all;
clear all;

% Set parameters
numDims = 40;
m = 10*numDims;
numDimsPerGroup = 4;

% Test 1
% 2 1D examples
fprintf('Test 1\n======================\n');
[f, fProps] = getAdditiveFunction(numDims, numDimsPerGroup);
decomposition = fProps.decomposition;
numGroups = numDims/numDimsPerGroup;
bounds = fProps.bounds;

% Generate Data
X = bsxfun(@plus, bsxfun(@times, rand(m, numDims), ...
          (bounds(:,2) - bounds(:,1))' ), bounds(:,1)');
Y = f(X);

% Perform Additive GP Regression
hyperParams.useSameSm = true; %false;
hyperParams.useSamePr = true;
hyperParams.fixSm = false;
hyperParams.fixPr = false;
hyperParams.sigmaSmRanges = [];
hyperParams.sigmaSmRange = [];
hyperParams.sigmaPrRanges = [];
hyperParams.sigmaPrRange = [];
hyperParams.noises = 0.00 * std(Y) * ones(numGroups, 1);
hyperParams.commonNoise = 0.01 * std(Y);
hyperParams.meanFuncs = [];
hyperParams.commonMeanFunc = [];
dummyPt = zeros(0, numDims);
[mu, KPost, Mus, KPosts, combinedFuncH, funcHs, sigmaSmOpts, sigmaPrOpts] = ...
  addGPMargLikelihood(X, Y, dummyPt, decomposition, hyperParams);
fprintf('Chosen Smoothness: %s\n', mat2str(sigmaSmOpts));
fprintf('Chosen Prior: %s\n', mat2str(sigmaPrOpts));
% % Now plot the functions
% plotAddGPOutputs(combinedFuncH, funcHs, f, trueFuncHs, decomposition, bounds);

% Now compare with a generic GP Regression
hp.meanFunc = [];
hp.sigmaSmRange = [];
hp.sigmaPrRange = [];
hp.sigmaSm = 0;
hp.sigmaPr = 0;
hp.noise = hyperParams.commonNoise;
hp.retFunc = true;
[~, ~, naiveGPFuncH] = GPMargLikelihood(X, Y, dummyPt, hp);

% Now do the test
numTestPts = 100*numDims;
testPts =  bsxfun(@plus, bsxfun(@times, rand(numTestPts, numDims), ...
               (bounds(:,2) - bounds(:,1))' ), bounds(:,1)');
trueVals = f(testPts);
addGPVals = combinedFuncH(testPts);
addGPErr = norm(trueVals-addGPVals)/sqrt(numTestPts); 
naiveGPVals = naiveGPFuncH(testPts);
naiveGPErr = norm(trueVals-naiveGPVals)/sqrt(numTestPts); 
fprintf('Add-GP Error: %0.5f\nNaive-GP Error:%0.5f\n', ...
  addGPErr, naiveGPErr);

