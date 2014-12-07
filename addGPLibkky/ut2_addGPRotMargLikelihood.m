% Unit Test for additive GP Regression

addpath ~/libs/kky-matlab/GPLibkky/
addpath ~/libs/kky-matlab/utils/
addpath ~/libs/kky-matlab/ancillary/
addpath ../BOLibkky/
addpath ../FOptM/
addpath ..//
addpath ../utils/

close all;
clear all;
clc; clc;

% Set parameters
numDims = 9;
numDimsPerGroup = 3;
m = 60*numDimsPerGroup;

% Test 1
% 2 1D examples
fprintf('Test 1\n======================\n');
[f, fProps] = getAdditiveFunction(numDims, numDimsPerGroup);
decomposition = fProps.decomposition;
numGroups = floor(numDims/numDimsPerGroup);
bounds = fProps.bounds;

% Generate Data
X = bsxfun(@plus, bsxfun(@times, rand(m, numDims), ...
          (bounds(:,2) - bounds(:,1))' ), bounds(:,1)');
Y = f(X);

% Perform Additive GP Regression with unknown decomposition
uhps.sigmaSmRange = [];
uhps.sigmaPrRange = [];
uhps.noises = 0;
uhps.commonNoise = 0.01 * std(Y);
uhps.meanFuncs = [];
uhps.commonMeanFunc = [];
uhps.numAOptInits = 1;
uhps.numAOptIters = 10;
uhps.numOuterInits = 1;
uhps.numBwSigmaDiRectIters = [];
uhps.decompStrategy = 'learn';
uhps.decompStrategy = 'partialLearn';
dummyPt = zeros(0,numDims);
[~, ~, ~, ~, combinedXFuncH, combinedZFuncH, funcHs, ...
  sigmaSmOpts, sigmaPrOpts, A, learnedDecomp] = ...
  addGPRotMargLikelihood(X, Y, dummyPt, numDimsPerGroup, numGroups, uhps);
fprintf('Chosen Smoothness: %s\n', mat2str(sigmaSmOpts));
fprintf('Chosen Prior: %s\n', mat2str(sigmaPrOpts));
if size(A, 1) > 10, A(1:10, 1:10),
else, A,
end
[~, permutes] = orthToPermutation(A),
fProps.permuteOrder,

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
ukAddGPVals = combinedXFuncH(testPts);
ukAddGPErr = norm(trueVals - ukAddGPVals)/sqrt(numTestPts);
addGPVals = combinedFuncH(testPts);
addGPErr = norm(trueVals-addGPVals)/sqrt(numTestPts); 
naiveGPVals = naiveGPFuncH(testPts);
naiveGPErr = norm(trueVals-naiveGPVals)/sqrt(numTestPts); 
fprintf('UK Add-GP Error: %0.5f\nAdd-GP Error: %0.5f\nNaive-GP Error:%0.5f\n', ...
  ukAddGPErr, addGPErr, naiveGPErr);

