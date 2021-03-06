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

trial = 2;
uTest = false;
% uTest = true;

switch trial
  case 1
    numDims = 10;
    trueNumDimsPerGroup = 3;
    numDimsPerGroupCands = [1 3 5 10]';
  case 2 
    numDims = 24;
    trueNumDimsPerGroup = 12;
    numDimsPerGroupCands = [1 4 6 12 24]';
  case 3
    numDims = 40;
    trueNumDimsPerGroup = 18;
    numDimsPerGroupCands = [1 4 10 20 40]';
  case 4
    numDims = 60;
    trueNumDimsPerGroup = 25;
    numDimsPerGroupCands = [1 6 10 25 60]';
  case 5
    numDims = 96;
    trueNumDimsPerGroup = 29;
    numDimsPerGroupCands = [4 8 16 32 96]';
  case 6
    numDims = 120;
    trueNumDimsPerGroup = 55;
    numDimsPerGroupCands = [8 15 30 55 120]';
  otherwise
end
numdCands = numel(numDimsPerGroupCands);

if ~uTest
  % Fixed experiment parameters
  % numExperiments = 3;
  numExperiments = 20;
  % numDiRectEvals = 500;
  % numDiRectEvals = min(5000, max(100*numDims, 500));
  numDiRectEvals = min(5000, max(40*numDims, 500));
  numIters = 830;
else
  % Unit text experiment parameters
  numExperiments = 1;
  numDiRectEvals = 5;
  numIters = 10;
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
saveFileName = sprintf('%stoyExp-%d-%d-%s.mat', resultsDir, numDims, ...
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
% use the same initialization points
boParams.initPts = rand(boParams.numInitPts, numDims);
boParams.initVals = func(boParams.initPts);
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
boParams.diRectParams.maxevals = ceil(numDiRectEvals/trueNumGroups);
boParams.diRectParams.maxits = inf;

totalNumQueries = numIters + boParams.numInitPts

% From here on customize each parameters separately.
% Initialize optimization parameters
boAddParams = boParams;
boAddParams.decompStrategy = 'partialLearn';

boUDParams = boParams;
boUDParams.decompStrategy = 'decide';
choosedM = {'mll', 'sample', 'rand'};
numChoosedM = numel(choosedM);

% Initialize arrays for storing the history
boAddHistories = zeros(numExperiments, totalNumQueries, numdCands);
boUDHistories = zeros(numExperiments, totalNumQueries, numChoosedM);

% For storing simple regret values
boAddSimpleRegrets = zeros(numExperiments, totalNumQueries, numdCands);
boUDSimpleRegrets = zeros(numExperiments, totalNumQueries, numChoosedM);
randSimpleRegrets = zeros(numExperiments, totalNumQueries);
% For storing cumulative regret values
boAddCumRegrets = zeros(numExperiments, totalNumQueries, numdCands);
boUDCumRegrets = zeros(numExperiments, totalNumQueries, numChoosedM);
randCumRegrets = zeros(numExperiments, totalNumQueries);
MHistAll = [];
ptAll = [];


% Initialize (d,M) pairs that passed to bayesOptDecideAddGP
decomp = cell(numdCands,1);
for i = 1:numdCands
  decomp{i}.d = numDimsPerGroupCands(i);
  decomp{i}.M = floor(numDims / decomp{i}.d);
end

for expIter = 1:numExperiments

  fprintf('\n==============================================================\n');
  fprintf('Experiment %d/ %d\nMaxVal: %0.4f\n', ...
    expIter, numExperiments, trueMaxVal);
  fprintf('Num DiRectEvals: %d\n', numDiRectEvals);
  fprintf('==============================================================\n');

  % Learn Decomposition
  fprintf('\nKnown Grouping Unknown Decomposition\n');

  % For the candidates in numDimsPerGroupCands
  for candIter = 1:numdCands
    fprintf('\nUsing an arbitrary %d/ %d decomposition\n', ...
      numDimsPerGroupCands(candIter), numDims );

    decompAdd = decomp{candIter};
    [~, ~, ~, valHistAdd] = ...
      decide(func, decompAdd, bounds, numIters, boAddParams);

    [sR, cR] = getRegrets(trueMaxVal, valHistAdd);
    boAddHistories(expIter, :, candIter) = valHistAdd';
    boAddSimpleRegrets(expIter, :, candIter) = sR';
    boAddCumRegrets(expIter, :, candIter) = cR';
  end


  % Decomposition varies across iterations  
  for iter=1:numChoosedM
    fprintf('\nChoose the decomposition\n');
    boUDParams.choosedM = choosedM{iter};

    [~, ~, ~, valHistDecide] = ...
      decide(func, decomp, bounds, numIters, boUDParams);

    [sR, cR] = getRegrets(trueMaxVal, valHistDecide);
    boUDHistories(expIter, :, iter) = valHistDecide';
    boUDSimpleRegrets(expIter, :, iter) = sR';
    boUDCumRegrets(expIter, :, iter) = cR';
  end

  fprintf('\nRandom optimization\n');
  randQueries = bsxfun(@plus, ...
    bsxfun(@times, rand(totalNumQueries, numDims), ...
    (bounds(:,2) - bounds(:,1))' ), bounds(:,1)' );
  randHistories(expIter, :) = func(randQueries)';
  [sR, cR] = getRegrets(trueMaxVal, randHistories(expIter, :)');
  randSimpleRegrets(expIter, :) = sR';
  randCumRegrets(expIter, :) = cR';


  % Save Results at each iteration
  if ~uTest
    save(saveFileName, 'numDims', 'trueNumDimsPerGroup', 'func', ...
      'funcProperties', 'trueMaxVal',  ...
      'numIters', 'totalNumQueries', ...
      'boAddHistories', 'boUDHistories', ...
      'boAddSimpleRegrets', 'boUDSimpleRegrets', ...
      'boAddCumRegrets', 'boUDCumRegrets', 'boUDParams', ...
      'randSimpleRegrets', 'randCumRegrets', ...
      'numDimsPerGroupCands', 'choosedM');
  end
end
