%Unit test for bayesOptGP.m

close all;
clear all;
addpath ~/libs/kky-matlab/utils/
addpath ~/libs/kky-matlab/ancillary/
addpath ~/libs/kky-matlab/GPLibkky/
addpath ../BOLibkky/
addpath ../addGPLibkky/
addpath ../utils/

warning off;

% Choose initial parameters
numExperiments = 3;
numDims = 20;
numDimsPerGroup = 4;
numGroups = floor(numDims/numDimsPerGroup);

% Obtain the Function
[func, funcProperties] = getAdditiveFunction(numDims, numDimsPerGroup);
bounds = funcProperties.bounds;
numIters = min(8 * 2^numDims, 1000);
maxPt = funcProperties.maxPt;
maxVal = funcProperties.maxVal;
decomposition = funcProperties.decomposition;
fprintf('True maxVal, maxPt: %f, %s\n', maxVal, mat2str(maxPt) );

% Set parameters for Regular Bayesian Optimization
boParams.optPtStdThreshold = 0.002;
boParams.alBWLB = 1e-5;
boParams.alBWUB = 5;
boParams.numInitPts = 4*numDims;
boParams.gpNoiseLevel = 0.5; % * std(f(th));
boParams.utilityFunc = 'UCB';

% Set parameters for Additive Bayesian Optimization
% First the common params
boAddParams.optPtStdThreshold = boParams.optPtStdThreshold;
boAddParams.alBWLB = boParams.alBWLB;
boAddParams.alBWUB = boParams.alBWUB;
boAddParams.numInitPts = boParams.numInitPts;
boAddParams.commonNoise = boParams.gpNoiseLevel;
boAddParams.utilityFunc = boParams.utilityFunc;
% Additional Params
boAddParams.meanFuncs = [];
boAddParams.commonMeanFunc = [];
boAddParams.noises = 0 * ones(numGroups, 1);
boAddParams.fixPr = false;
boAddParams.useSamePr = true;
boAddParams.sigmaPrRange = []; % Let the Marginal Likelihood pick this
boAddParams.useFixedBandWidth = false;
boAddParams.fixSm = false;
boAddParams.useSameSm = true;

% Set parameters for individual additive Bayesian Optimization
boIndParams.optPtStdThreshold = boParams.optPtStdThreshold;
boIndParams.alBWLB = boParams.alBWLB;
boIndParams.alBWUB = boParams.alBWUB;
boIndParams.numInitPtsPerGroup = round(boParams.numInitPts/numGroups);
boIndParams.noise = boParams.gpNoiseLevel;
boIndParams.utilityFunc = boParams.utilityFunc;
boIndParams.sigmaPrRange = [];
boIndParams.meanFunc = [];
% aDditional params
% boIndParams.meanFuncs = [];
% boIndParams.commonMeanFunc = [];
% boIndParams.noises = 0 * ones(numGroups, 1);
% boIndParams.fixPr = false;
% boIndParams.useSamePr = true;
% boIndParams.sigmaPrRange = []; % Let the Marginal Likelihood pick this
% boIndParams.useFixedBandWidth = false;
% boIndParams.fixSm = false;
% boIndParams.useSameSm = true;

% Initialize arrays for storing the history
totalNumQueries = numIters + boParams.numInitPts;
boHistories = zeros(numExperiments, totalNumQueries);
boIndHistories = zeros(numExperiments, totalNumQueries);
boAddHistories = zeros(numExperiments, totalNumQueries);
randHistories = zeros(numExperiments, totalNumQueries);

for expIter = 1:numExperiments

  % Call Additive Individual BO
  [~, ~, ~, ~, boIndHist] = ...
    bayesOptIndAddGP(func, decomposition, bounds, numIters, boIndParams);
  boIndHistories(expIter, :) = boIndHist';

  % Call Regular BO
  [boMaxVal, boMaxPt, boQueries, boVals, boHist] = ...
    bayesOptGP(func, bounds, numIters, boParams);
  boHistories(expIter, :) = boHist';

  % Call Additive BO
  [boAddMaxVal, boAddMaxPt, boAddQueries, boAddVals, boAddHist] = ...
    bayesOptAddGP(func, decomposition, bounds, numIters, boAddParams);
  boAddHistories(expIter, :) = boAddHist';

  % Random calls
  randQueries = bsxfun(@plus, ...
          bsxfun(@times, rand(totalNumQueries, numDims), ...
                 (bounds(:,2) - bounds(:,1))' ), bounds(:,1)');
  randVals = func(randQueries);
  randHistories(expIter, :) = max(randVals(cumsum(triu(ones(length(randVals))))));

end

  % finally Call DIRECT
  diRectOpts.maxevals = totalNumQueries;
  [~, ~, diRectHist] = diRectWrap(func, bounds, diRectOpts);
  numFuncEvaluations = diRectHist(end, 2);

  boHistMean = mean(boHistories, 1); boHistStd = std(boHistories, 1);
  boAddHistMean = mean(boAddHistories, 1); boAddHistStd = std(boAddHistories, 1);
  boIndHistMean = mean(boIndHistories, 1); boIndHistStd = std(boIndHistories, 1);
  randHistMean = mean(randHistories, 1); randHistStd = std(randHistories, 1);

  % 4. Plot iteration statistics
  figure;  hold on,
  loglog(1:totalNumQueries, boHistMean, 'b'); hold on,
  loglog(1:totalNumQueries, boAddHistMean, 'c'); hold on;
  loglog(1:totalNumQueries, boIndHistMean, 'g')
  loglog(diRectHist(:,2),diRectHist(:,3), 'r');
  loglog(1:totalNumQueries, randHistMean, 'k');
  legend('BO', 'AddBO', 'AddIndBO', 'DiRect', 'Rand');
  if numExperiments > 1
  errorbar(1:totalNumQueries, boHistMean, boHistStd/sqrt(numExperiments), 'Color', 'b');
  errorbar(1:totalNumQueries, boAddHistMean, boAddHistStd/sqrt(numExperiments), ...
          'Color', 'c');
  errorbar(1:totalNumQueries, boIndHistMean, boIndHistStd/sqrt(numExperiments), ...
          'Color', 'g');
  errorbar(1:totalNumQueries, randHistMean, randHistStd/sqrt(numExperiments), 'Color', 'k');
  end
  loglog([0, totalNumQueries], [maxVal, maxVal], 'k')';
%   axis([20*numDims totalNumQueries*1.1, maxVal - (maxVal-randHistMean(1))/1, maxVal*1.1]);

