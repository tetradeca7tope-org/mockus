%BO for LRGs

close all;
clear all;
LIBPATH = '../add-gp-bandits-v1/';
addpath(genpath(LIBPATH));
warning off;

% Problem parameters
numExperiments = 2;
numIters = 100;
numDims = 15;  numDimsPerGroupCands = [15;1;2;3;5];

numInitPts = 10;
numDiRectEvals = 2000; 

% Other derived Parameters
numdCands = numel(numDimsPerGroupCands);
bounds = repmat([0 1], numDims, 1);

% Get the function
eapExp = EAPExperiment();
func = @(arg) eapExp.normCoordFitness(arg);

% Ancillary stuff
resultsDir = 'results/';
saveFileName = sprintf('%seap%d-%s-%s.mat', resultsDir, numDims,...
  mat2str(numDimsPerGroupCands), datestr(now,'ddmm-hhMMss') );
saveFileName,

% Parameters for additive Bayesian optimisation
boParams.optPtStdThreshold = 0.002;
boParams.alBWLB = 1e-5;
boParams.alBWUB = 5;
boParams.commonNoise = 0.5;
boParams.utilityFunc = 'UCB';
boParams.meanFuncs = [];
boParams.commonMeanFunc = @(arg) zeros(size(arg, 1), 1);
% boParams.commonMeanFunc = []; %@(arg) meanFth * ones(size(arg, 1), 1);
boParams.useSamePr = true;
boParams.useSameSm = true;
boParams.fixPr = false;
boParams.fixSm = false;
boParams.sigmaPrRange = [0.03 30] * stdFth;
boParams.useFixedBandwidth = false;

% The rest - arbitrary decompositions
boAddParams = boParams;
boAddParams.decompStrategy = 'partialLearn';
boAddParams.diRectParams.maxits = inf;

totalNumQueries = numIters + numInitPts;
% Initialize an array for storing the query points
boAddQueryPts = zeros(totalNumQueries, numDims, numExperiments, numdCands);
randQueryPts = zeros(totalNumQueries, numDims, numExperiments);
% Initialize arrays for storing the histories
boAddHistories = zeros(numExperiments, totalNumQueries, numdCands);
randHistories = zeros(numExperiments, totalNumQueries);
% For storing maximum values
boAddMaxVals = zeros(numExperiments, totalNumQueries, numdCands);
randMaxVals = zeros(numExperiments, totalNumQueries);
% For storing cumulative rewards
boAddCumRewards = zeros(numExperiments, totalNumQueries, numdCands);
randCumRewards = zeros(numExperiments, totalNumQueries);
% Store times
boAddTimes = zeros(numExperiments, numdCands);
randTimes = zeros(numExperiments);

% First call at the nominal value
nominalVal = func(nominalPt);

% % First Call Direct
diRectHist = [(1:totalNumQueries)' (1:totalNumQueries)' zeros(totalNumQueries,1)];
diRectHistories = 0;
diRectOptPt = 0;
diRectTime = 0;
  fprintf('First Running DiRect\n============================================\n');
  diRectOpts.maxevals = totalNumQueries;
  diRectOpts.maxits = inf;
  diRectOpts.showits = true;
  tic;
[~, ~, diRectHist, diRectQueries, diRectHistory] = ...
  diRectWrap(func, bounds, diRectOpts);
  diRectTime = toc;
[diRectMaxVals, diRectCumRewards] = getRegrets(-inf, diRectHistory);
fprintf('DiRectOpt = %.4f\n', diRectMaxVals(end));


for expIter = 1:numExperiments

  % First do Random
  randQueries = bsxfun(@plus, ...
    bsxfun(@times, rand(totalNumQueries, numDims), ...
      (bounds(:,2) - bounds(:,1))' ), bounds(:,1)' );
  tic;
  randHistories(expIter, :) = func(randQueries)';
  randQueryPts(:,:,expIter) = randQueries;
  randTimes(expIter) = toc;
  [sR, cR] = getRegrets(-inf, randHistories(expIter, :)');
  randMaxVals(expIter, :) = sR';
  randCumRewards(expIter, :) = cR';
  fprintf('Random Max: %0.5f\n', randMaxVals(expIter, totalNumQueries));

  fprintf('\n\nFirst obtaining initialization for BO\n');
%   boAddParams.initPts = [0.1*rand(1, numDims) + nominalPt; ...
%     rand(numInitPts-1, numDims)];
  boAddParams.initPts = [rand(numInitPts, numDims)];
%   boAddParams.initPts = diRectQueries(1:numInitPts, :) + ...
%     0.01*rand(numInitPts, numDims);
  boAddParams.initVals = func(boAddParams.initPts);
  fprintf('Max Init Val: %0.5f\n', max(boAddParams.initVals));

  fprintf('\n==============================================================\n');
  fprintf('Experiment %d/ %d\nNominal Value: %.4f\n', ...
    expIter, numExperiments, nominalVal);
  fprintf('Num DiRectEvals: %d\n', numDiRectEvals);
  fprintf('==============================================================\n');

  % For the candidates in numDimsPerGroupCands
  for candIter = 1:numdCands
    fprintf('\nUsing a %d/ %d decomposition\n', ...
      numDimsPerGroupCands(candIter), numDims );
    [decompAdd, boAddParamsCurr, numCurrGroups] = ...
      getDecompForParams(numDims, numDimsPerGroupCands(candIter), ...
                  boAddParams, true);
    boAddParamsCurr.diRectParams.maxevals = ...
      ceil(numDiRectEvals/(numDims/numDimsPerGroupCands(candIter)));
    tic;
    [~, ~, queries, valHistAdd] = ...
      bayesOptDecompAddGP(func, decompAdd, bounds, numIters, boAddParamsCurr);
    boAddTimes(expIter, candIter) = toc;
    [sR, cR] = getRegrets(-inf, valHistAdd);
    boAddQueryPts(:,:,expIter, candIter) = queries;
    boAddHistories(expIter, :, candIter) = valHistAdd';
    boAddMaxVals(expIter, :, candIter) = sR';
    boAddCumRewards(expIter, :, candIter) = cR';
  end

  % Save Results at each iteration
  save(saveFileName, 'numDims', 'numDimsPerGroupCands', ...
    'numIters', 'totalNumQueries', 'numdCands', ...
    'boAddHistories', 'randHistories', 'diRectHistory', ...
    'boAddMaxVals', 'randMaxVals', 'diRectMaxVals', ...
    'boAddCumRewards', 'randCumRewards', 'diRectCumRewards' , ...
    'boAddQueryPts', 'randQueryPts', 'diRectOptPt', ...
    'boAddTimes', 'randTimes', 'diRectTime');
  
end

  beep;
%   clearvars -except saveFileName;
  load(saveFileName);
  plotLRGResults;
