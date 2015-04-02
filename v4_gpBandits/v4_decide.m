% Experiment Set up for Bayesian Optimization and GP Bandits

close all;
clear all;
LIBPATH = '~/libs/kky-matlab/';
addpath([LIBPATH 'utils']); 
addpath([LIBPATH 'ancillary']); 
addpath([LIBPATH 'GPLibkky']); 
addpath ../BOLibkky/
addpath ../addGPLibkky/
addpath ../utils/

warning off;

% Problem parameters
numExperiments = 2;
numDims = 40;

numDiRectEvals = 500;
trueNumDimsPerGroup = 8;
numIters = 400;

numDimsPerGroupCands = [40 4 8 10 16]';

% tiny experiments 
%numExperiments = 1;
%numDims = 4;
%
%numDiRectEvals = 5;
%trueNumDimsPerGroup = 2;
%numIters = 12;
%
%numDimsPerGroupCands = [4 1 2]';

numdCand = numel(numDimsPerGroupCands);

% Initialize (d,M) pairs that passed to bayesOptDecideAddGP
decomp = cell(numdCands-1,1);
for i = 1:numdCands-1
  decomp{i}.d = numDimsPerGroupCands(i+1);
  decomp{i}.M = floor(numDims / decomp{i}.d);
end

% Get the function
[func, funcProperties] = getAdditiveFunction(numDims, trueNumDimsPerGroup);
bounds = funcProperties.bounds;
trueDecomp = funcProperties.decomposition;
trueMaxVal = funcProperties.maxVal;
trueMaxPt = funcProperties.maxPt;
trueNumGroups = numel(trueDecomp);

% Ancillary stuff
resultsDir = 'results/';
saveFileName = sprintf('%s-toyExp-%d-%d-%s.mat', resultsDir, numDims, ...
  trueNumDimsPerGroup, datestr(now,'ddmm-hhMMss') );

% Compute some statistics to help with the initialization
th = rand(1000, numDims); fth = func(th);
meanFth = mean(fth);
stdFth = std(fth);

% Parameters for additive Bayesian optimization
boParams.optPtStdThreshold = 0.002;
boParams.alBWLB = 1e-5;
boParams.alBWUB = 5;
boParams.numInitPts = 10; 
boParams.commonNoise = 0.01 * stdFth;
boParams.utilityFunc = 'UCB';
boParams.meanFuncs = [];
boParams.commonMeanFunc = @(arg) zeros(size(arg, 1), 1);
boParams.useSamePr = true;
boParams.useSameSm = true;
boParams.fixPr = false;
boParams.fixSm = false;
boParams.sigmaPrRange = [0.03 30] * stdFth;
boParams.useFixedBandwidth = false;


totalNumQueries = numIters + boParams.numInitPts


% From here on customize each parameters separately.
% Initialize optimization parameters
boKDParams = boParams;
boKDParams.decompStrategy = 'known';
boKDParams.diRectParams.maxevals = ceil(numDiRectEvals/trueNumGroups);
boKDParams.diRectParams.maxits = inf;

boAddParams = boParams;
boAddParams.decompStrategy = 'partialLearn';

boUDParams = boKDParams;
boUDParams.decompStrategy = 'decide';

% Initialize arrays for storing the history
boKDHistories = zeros(numExperiments, totalNumQueries);
boAddHistories = zeros(numExperiments, totalNumQueries, numdCands);
boUDHistories = zeros(numExperiments, totalNumQueries);

% For storing simple regret values
boKDSimpleRegrets = zeros(numExperiments, totalNumQueries);
boAddSimpleRegrets = zeros(numExperiments, totalNumQueries, numdCands);
boUDSimpleRegrets = zeros(numExperiments, totalNumQueries);

% For storing cumulative regret values
boKDCumRegrets = zeros(numExperiments, totalNumQueries);
boAddCumRegrets = zeros(numExperiments, totalNumQueries, numdCands);
boUDCumRegrets = zeros(numExperiments, totalNumQueries);


for expIter = 1:numExperiments

  fprintf('\n==============================================================\n');
  fprintf('Experiment %d/ %d\nMaxVal: %0.4f\n', ...
    expIter, numExperiments, trueMaxVal);
  fprintf('Num DiRectEvals: %d\n', numDiRectEvals);
  fprintf('==============================================================\n');

  % Known true decomposition
  fprintf('\nKnown Decomposition\n');
  boKDParams.noises = 0 * ones(trueNumGroups, 1);
  [~, ~, ~, valHistKD] = ...
    bayesOptDecompAddGP(func, trueDecomp, bounds, numIters, boKDParams);
  [sR, cR] = getRegrets(trueMaxVal, valHistKD);

  boKDHistories(expIter, :) = valHistKD';
  boKDSimpleRegrets(expIter, :) = sR';
  boKDCumRegrets(expIter, :) = cR';

  % Learn Decomposition
  fprintf('\nKnown Grouping Unknown Decomposition\n');

  % For the candidates in numDimsPerGroupCands
  for candIter = 1:numdCands
    fprintf('\nUsing an arbitrary %d/ %d decomposition\n', ...
      numDimsPerGroupCands(candIter), numDims );

    [decompAdd, boAddParamsCurr, numCurrGroups] = ...
      getDecompForParams(numDims, numDimsPerGroupCands(candIter), ...
                  boAddParams, true);
    
    boAddParamsCurr.diRectParams.maxevals = ceil(0.9 * numDiRectEvals/numCurrGroups);
    
    [~, ~, ~, valHistAdd] = ...
      bayesOptDecompAddGP(func, decompAdd, bounds, numIters, boAddParamsCurr);
    
    fprintf('In additive model, size of valHistAdd %d\n', numel(valHistAdd));

    [sR, cR] = getRegrets(trueMaxVal, valHistAdd);
    boAddHistories(expIter, :, candIter) = valHistAdd';
    boAddSimpleRegrets(expIter, :, candIter) = sR';
    boAddCumRegrets(expIter, :, candIter) = cR';
 
  end

  
  % Decomposition varies across iterations
  fprintf('\nDecomposition varies across iterations\n');

  [~, ~, ~, valHistDecide] = ...
      bayesOptDecideAddGP(func, decomp, bounds, numIters, boUDParams); 
  
  [sR, cR] = getRegrets(trueMaxVal, valHistDecide);
  
  size(valHistDecide);
  size(boUDHistories);
  boUDHistories(expIter, :) = valHistDecide';
  boUDSimpleRegrets(expIter, :) = sR';
  boUDCumRegrets(expIter, :) = cR';

  % Save Results at each iteration
  save(saveFileName, 'numDims', 'trueNumDimsPerGroup', 'func', ...
    'funcProperties', 'trueMaxVal',  ...
    'numIters', 'totalNumQueries', ...
    'boKDHistories', 'boAddHistories', 'boUDHistories', ...
    'boKDSimpleRegrets', 'boAddSimpleRegrets', 'boUDSimpleRegrets', ...
    'boKDCumRegrets', 'boAddCumRegrets', 'boUDCumRegrets', ...
    'numDimsPerGroupCands');

end
