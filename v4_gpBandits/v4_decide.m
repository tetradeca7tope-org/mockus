% Experiment Set up for Bayesian Optimization and GP Bandits

close all;
clear all;
LIBPATH = '../../libs/';
addpath([LIBPATH 'utils']); 
addpath([LIBPATH 'ancillary']); 
addpath([LIBPATH 'GPLibkky']); 
addpath ../BOLibkky/
addpath ../addGPLibkky/
addpath ../utils/

warning off;

% Problem parameters
numExperiments = 1;
numDims = 10;

numDiRectEvals = 500;
trueNumDimsPerGroup = 3;
numIters = 100;


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

% From here on customize each parameters separately.
% Known Decomposition
boKDParams = boParams;
boKDParams.decompStrategy = 'known';
boKDParams.diRectParams.maxevals = ceil(numDiRectEvals/trueNumGroups);
boKDParams.diRectParams.maxits = inf;

% Known grouping but not decomposition
boUDParams = boKDParams;
boUDParams.decompStrategy = 'decide';

decomp = cell(4,1);
decomp{1}.d=1;
decomp{1}.M=10;

decomp{2}.d=2;
decomp{2}.M=5;

decomp{3}.d=5;
decomp{3}.M=2;

decomp{4}.d=10;
decomp{4}.M=1;


totalNumQueries = numIters + boParams.numInitPts;

% Initialize arrays for storing the history
boKDHistories = zeros(numExperiments, totalNumQueries);
boUDHistories = zeros(numExperiments, totalNumQueries);

% For storing simple regret values
boKDSimpleRegrets = zeros(numExperiments, totalNumQueries);
boUDSimpleRegrets = zeros(numExperiments, totalNumQueries);

% For storing cumulative regret values
boKDCumRegrets = zeros(numExperiments, totalNumQueries);
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

  % Learn the decomposition 
  fprintf('\nLearn Decomposition\n');
  [~, ~, ~, valHistAdd] = ...
      bayesOptDecideAddGP(func, decomp, bounds, numIters, boUDParams); 
  
  [sR, cR] = getRegrets(trueMaxVal, valHistAdd);
  
  % TODO: dimensions dis-match bug:
  size(valHistAdd');
  size(boUDHistories);

 % boUDHistories(expIter, :) = valHistAdd';
 % boUDSimpleRegrets(expIter, :) = sR';
 % boUDCumRegrets(expIter, :) = cR';

  % Save Results at each iteration
%   save(saveFileName, 'numDims', 'trueNumDimsPerGroup', 'func', ...
%    'funcProperties', 'trueMaxVal',  ...
%    'numIters', 'totalNumQueries', ...
%    'boKDHistories', 'boUDHistories', ...
%    'boKDSimpleRegrets', 'boUDSimpleRegrets', ...
%    'boKDCumRegrets', 'boUDCumRegrets');

end
