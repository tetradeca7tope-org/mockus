% testing the BO stuff on non-additive functions

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
numExperiments = 4;
numDims = 24;
% numDims = 4;
% numDimsPerGroup = 2;

% Obtain the Function
[func, funcProperties] = getAdditiveFunction(numDims, 7);
funcProperties,
funcProperties.decomposition{1},
bounds = funcProperties.bounds;
numIters = min(8 * 2^numDims, 500);
maxPt = funcProperties.maxPt;
maxVal = funcProperties.maxVal;
fprintf('True maxVal, maxPt: %f, %s\n', maxVal, mat2str(maxPt) );

% Set parameters for regular Bayesian Optimization
boParams.optPtStdThreshold = 0.002;
boParams.alBWLB = 1e-5;
boParams.alBWUB = 5;
boParams.numInitPts = min(50, numDims);
boParams.gpNoiseLevel = 0.5; % * std(f(th));
boParams.utilityFunc = 'UCB';

% First the common params
boAddParams.optPtStdThreshold = boParams.optPtStdThreshold;
boAddParams.alBWLB = boParams.alBWLB;
boAddParams.alBWUB = boParams.alBWUB;
boAddParams.numInitPts = boParams.numInitPts;
boAddParams.commonNoise = boParams.gpNoiseLevel;
boAddParams.utilityFunc = 'UCB'; %boParams.utilityFunc;
% Additional Params
boAddParams.meanFuncs = [];
boAddParams.commonMeanFunc = @(arg) zeros(size(arg,1), 1); 
boAddParams.fixPr = false;
boAddParams.useSamePr = true;
boAddParams.sigmaPrRange = []; % Let the Marginal Likelihood pick this
boAddParams.useFixedBandWidth = false;
boAddParams.fixSm = false;
boAddParams.useSameSm = true;
% boAddParams.knowDecomp = true;
boAddParams.decompStrategy = 'partialLearn';

% Initialize arrays for storing the history
totalNumQueries = numIters + boParams.numInitPts;
boHistories = zeros(numExperiments, totalNumQueries);
boAddHistories1 = zeros(numExperiments, totalNumQueries);
boAddHistories2 = zeros(numExperiments, totalNumQueries);
boAddHistories3 = zeros(numExperiments, totalNumQueries);
randHistories = zeros(numExperiments, totalNumQueries);

% numDimsPerGroupVals = [1; 5; 10];
numDimsPerGroupVals = [1; 6; 12];

for expIter = 1:numExperiments

  fprintf('Experiment %d/ %d\n', expIter, numExperiments);
  fprintf('==================================================================');

  % Call Regular BO
  fprintf('\nRegular BO\n');
  numDimsPerGroup = numDims;
  numGroups = 1; 
  boParams = boAddParams;
  boParams.decompStrategy = 'known';
  boAddParams.noises = 0;
  decomp0 = cell(1,1); decomp0{1} = 1:numDims;
  [boMaxVal, boMaxPt, boQueries, boVals, boHist] = ...
    bayesOptUDAddGP(func, decomp0, bounds, numIters, boParams);
  boHistories(expIter, :) = boHist';

  % Call Additive BO with different numDims per Group 
  fprintf('\nAdditive BO 1\n');
  numDimsPerGroup = numDimsPerGroupVals(1);
  numGroups = floor(numDims/numDimsPerGroup);
  if strcmp(boAddParams.decompStrategy, 'known')
    boAddParams.noises = 0 * ones(numGroups, 1);
    decomp1 = cell(numGroups, 1);
    for i = 1:numGroups
      decomp1{i} = ( (i-1)*numDimsPerGroup + 1 : i*numDimsPerGroup );
    end
  else
    decomp1.d = numDimsPerGroup;
    decomp1.M = numGroups;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [boAddMaxVal, boAddMaxPt, boAddQueries, boAddVals, boAddHist] = ...
    bayesOptUDAddGP(func, decomp1, bounds, numIters, boAddParams);
  boAddHistories1(expIter, :) = boAddHist';

  % Call Additive BO  with different numDims per Group 
  fprintf('\nAdditive BO 2\n');
  numDimsPerGroup = numDimsPerGroupVals(2);
  numGroups = floor(numDims/numDimsPerGroup);
  if strcmp(boAddParams.decompStrategy, 'known')
    boAddParams.noises = 0 * ones(numGroups, 1);
    decomp2 = cell(numGroups, 1);
    for i = 1:numGroups
      decomp2{i} = ( (i-1)*numDimsPerGroup + 1 : i*numDimsPerGroup );
    end
  else
    decomp2.d = numDimsPerGroup;
    decomp2.M = numGroups;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [boAddMaxVal, boAddMaxPt, boAddQueries, boAddVals, boAddHist] = ...
    bayesOptUDAddGP(func, decomp2, bounds, numIters, boAddParams);
  boAddHistories2(expIter, :) = boAddHist';

  % Call Additive BO  with different numDims per Group 
  fprintf('\nAdditive BO 3\n');
  numDimsPerGroup = numDimsPerGroupVals(3);
  numGroups = floor(numDims/numDimsPerGroup);
  if strcmp(boAddParams.decompStrategy, 'known')
    boAddParams.noises = 0 * ones(numGroups, 1);
    decomp3 = cell(numGroups, 1);
    for i = 1:numGroups
      decomp3{i} = ( (i-1)*numDimsPerGroup + 1 : i*numDimsPerGroup );
    end
  else
    decomp3.d = numDimsPerGroup;
    decomp3.M = numGroups;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [boAddMaxVal, boAddMaxPt, boAddQueries, boAddVals, boAddHist] = ...
    bayesOptUDAddGP(func, decomp3, bounds, numIters, boAddParams);
  boAddHistories3(expIter, :) = boAddHist';


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
  boAddHistMean1 = mean(boAddHistories1, 1); boAddHistStd1 = std(boAddHistories1, 1);
  boAddHistMean2 = mean(boAddHistories2, 1); boAddHistStd2 = std(boAddHistories2, 1);
  boAddHistMean3 = mean(boAddHistories3, 1); boAddHistStd3 = std(boAddHistories3, 1);
  randHistMean = mean(randHistories, 1); randHistStd = std(randHistories, 1);

  % 4. Plot iteration statistics
  figure;  hold on,
  loglog(1:totalNumQueries, boHistMean, 'b'); hold on,
  loglog(1:totalNumQueries, boAddHistMean1, 'c'); hold on;
  loglog(1:totalNumQueries, boAddHistMean2, 'k'); hold on;
  loglog(1:totalNumQueries, boAddHistMean3, 'g'); hold on;
  loglog(diRectHist(:,2), diRectHist(:,3), 'r');
  loglog(1:totalNumQueries, randHistMean, 'm');
    % Create the legend
    boLeg = sprintf('BO %0.4f', max(boHistMean));
    boAddLeg1 = sprintf('AddBO-%d %0.4f', numDimsPerGroupVals(1), max(boAddHistMean1));
    boAddLeg2 = sprintf('AddBO-%d %0.4f', numDimsPerGroupVals(2), max(boAddHistMean2));
    boAddLeg3 = sprintf('AddBO-%d %0.4f', numDimsPerGroupVals(3), max(boAddHistMean3));
    diRectLeg = sprintf('DiRect %0.4f', max(diRectHist(:,3)));
    randLeg = sprintf('Rand %0.4f', max(randHistMean));
  legend(boLeg, boAddLeg1, boAddLeg2, boAddLeg3, diRectLeg, randLeg);

  if numExperiments > 1
  errorbar(1:totalNumQueries, boHistMean, boHistStd/sqrt(numExperiments), 'Color', 'b');
  errorbar(1:totalNumQueries, boAddHistMean1, boAddHistStd1/sqrt(numExperiments), ...
          'Color', 'c');
  errorbar(1:totalNumQueries, boAddHistMean2, boAddHistStd2/sqrt(numExperiments), ...
          'Color', 'k');
  errorbar(1:totalNumQueries, boAddHistMean3, boAddHistStd3/sqrt(numExperiments), ...
          'Color', 'g');
  errorbar(1:totalNumQueries, randHistMean, randHistStd/sqrt(numExperiments), 'Color', 'm');
  end
  loglog([0, totalNumQueries], [maxVal, maxVal], 'k')';
%   axis([20*numDims totalNumQueries*1.1, maxVal - (maxVal-randHistMean(1))/1, maxVal*1.1]);
  titleStr = sprintf('D=%d, maxVal=%0.4f\n', numDims, funcProperties.maxVal);
  title(titleStr);


