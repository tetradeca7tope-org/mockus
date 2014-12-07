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
numExperiments = 2;
numDims = 10;
numDimsPerGroup = 3;
% numDims = 4;
% numDimsPerGroup = 2;
numGroups = floor(numDims/numDimsPerGroup);

% Obtain the Function
[func, funcProperties] = getAdditiveFunction(numDims, numDimsPerGroup);
bounds = funcProperties.bounds;
numIters = min(8 * 2^numDims, 500);
maxPt = funcProperties.maxPt;
maxVal = funcProperties.maxVal;
decomposition = funcProperties.decomposition;
fprintf('True maxVal, maxPt: %f, %s\n', maxVal, mat2str(maxPt) );

% Set parameters for Regular Bayesian Optimization
boParams.optPtStdThreshold = 0.002;
boParams.alBWLB = 1e-5;
boParams.alBWUB = 5;
boParams.numInitPts = min(20, numDims);
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
boAddParams.commonMeanFunc = @(arg) zeros(size(arg,1), 1);
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
% boIndParams.meanFunc = []; TODO: set to be the zero function here.
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
boAddUDHistories = zeros(numExperiments, totalNumQueries);
randHistories = zeros(numExperiments, totalNumQueries);

for expIter = 1:numExperiments

  fprintf('Experiment %d/ %d\n', expIter, numExperiments);
  fprintf('==================================================================');

%   % Call Additive Individual BO
%   fprintf('\nAdditive Individual BO\n');
%   [~, ~, ~, ~, boIndHist] = ...
%     bayesOptIndAddGP(func, decomposition, bounds, numIters, boIndParams);
%   boIndHistories(expIter, :) = boIndHist';

  % Call Regular BO
  fprintf('\nRegular BO\n');
%   [boMaxVal, boMaxPt, boQueries, boVals, boHist] = ...
%     bayesOptGP(func, bounds, numIters, boParams);
%   boHistories(expIter, :) = boHist';
  numDimsPerGroup = numDims;
  numGroups = 1; 
  boParams = boAddParams;
  boParams.decompStrategy = 'known';
  boAddParams.noises = 0;
  decomp0 = cell(1,1); decomp0{1} = 1:numDims;
  [boMaxVal, boMaxPt, boQueries, boVals, boHist] = ...
    bayesOptUDAddGP(func, decomp0, bounds, numIters, boParams);
  boHistories(expIter, :) = boHist';

  % Call Additive BO
  fprintf('\nAdditive BO\n');
  boAddParams.decompStrategy = 'known';
  [boAddMaxVal, boAddMaxPt, boAddQueries, boAddVals, boAddHist] = ...
    bayesOptUDAddGP(func, decomposition, bounds, numIters, boAddParams);
  boAddHistories(expIter, :) = boAddHist';

  % Call Additive but with unknown decomposition
  fprintf('\nAdditive UD BO\n');
  decompParams.d = numel(decomposition{1});
  decompParams.M = numel(decomposition);
  boAddUDParams = boAddParams;
  boAddUDParams.decompStrategy = 'learn';
  [boAddUDMaxVal, boAddUDMaxPt, boAddUDQueries, boAddUDVals, boAddUDHist] = ...
    bayesOptUDAddGP(func, decompParams, bounds, numIters, boAddUDParams);
  boAddUDHistories(expIter, :) = boAddUDHist';

  % Random calls
  fprintf('\nRAND\n');
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
  boAddUDHistMean = mean(boAddUDHistories, 1); boAddUDHistStd = std(boAddUDHistories, 1);

  % 4. Plot iteration statistics
  figure;  hold on,
  loglog(1:totalNumQueries, boHistMean, 'b'); hold on,
  loglog(1:totalNumQueries, boAddHistMean, 'c'); hold on;
  loglog(1:totalNumQueries, boAddUDHistMean, 'k');
  loglog(1:totalNumQueries, boIndHistMean, 'g')
  loglog(diRectHist(:,2),diRectHist(:,3), 'r');
  loglog(1:totalNumQueries, randHistMean, 'm');
    % Create the legend
    boLeg = sprintf('BO %0.4f', max(boHistMean));
    boAddLeg = sprintf('AddBO %0.4f', max(boAddHistMean));
    boAddUDLeg = sprintf('AddUDBO %0.4f', max(boAddUDHistMean));
    boAddIndLeg = sprintf('AddIndBO %0.4f', max(boIndHistMean));
    diRectLeg = sprintf('DiRect %0.4f', max(diRectHist(:,3)));
    randLeg = sprintf('Rand %0.4f', max(randHistMean));
  legend(boLeg, boAddLeg, boAddUDLeg, boAddIndLeg, diRectLeg, randLeg);

  if numExperiments > 1
  errorbar(1:totalNumQueries, boHistMean, boHistStd/sqrt(numExperiments), 'Color', 'b');
  errorbar(1:totalNumQueries, boAddHistMean, boAddHistStd/sqrt(numExperiments), ...
          'Color', 'c');
  errorbar(1:totalNumQueries, boAddUDHistMean, boAddUDHistStd/sqrt(numExperiments), ...
          'Color', 'k');
  errorbar(1:totalNumQueries, boIndHistMean, boIndHistStd/sqrt(numExperiments), ...
          'Color', 'g');
  errorbar(1:totalNumQueries, randHistMean, randHistStd/sqrt(numExperiments), 'Color', 'm');
  end
  loglog([0, totalNumQueries], [maxVal, maxVal], 'k')';
%   axis([20*numDims totalNumQueries*1.1, maxVal - (maxVal-randHistMean(1))/1, maxVal*1.1]);
  titleStr = sprintf('D=%d, maxVal=%0.4f\n', numDims, funcProperties.maxVal);
  title(titleStr);

