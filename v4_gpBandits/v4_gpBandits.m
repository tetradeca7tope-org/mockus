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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are the experiments we will run
% numDimsPerGroupCands will contain arbitrary initial random partitiions. Of 
% them, one of them will contain all numDims dimensions so this is eqt to just 
% running naive BO.
% Then we will also run one experiment which knows the true decomposition and
% another which knows the grouping sizes and attempts to find them.

% Problem parameters
numExperiments = 1;
numDims = 40;
% numDimsPerGroupCands = [10 1 2 3 4]';
% numDimsPerGroupCands = [4]';
% numDimsPerGroupCands = [10 1 2 4]';
% numDimsPerGroupCands = [24 1 3 6 12]';
numDimsPerGroupCands = [40 1 5 10 20]';
% numDimsPerGroupCands = [60 1 5 10 20]';
% numDimsPerGroupCands = [96 1 4 12 24]';
% numDimsPerGroupCands = [200 1 5 20 40]';
% numDimsPerGroupCands = [300 1 5 30 60]';
% numDimsPerGroupCands = [500 1 5 25 50]';

trueNumDimsPerGroup = 18;
% Experiment parameters
numIters = 300;
numDiRectEvals = min(5000, max(100*numDims, 500));

numdCands = numel(numDimsPerGroupCands);
% Get the function
[func, funcProperties] = getAdditiveFunction(numDims, trueNumDimsPerGroup);
bounds = funcProperties.bounds;
trueDecomp = funcProperties.decomposition;
trueMaxVal = funcProperties.maxVal;
trueMaxPt = funcProperties.maxPt;
trueNumGroups = numel(trueDecomp);

% Ancillary stuff
resultsDir = 'results/';
saveFileName = sprintf('%sgpb%d-%d-%s-%s.mat', resultsDir, numDims, ...
  trueNumDimsPerGroup, mat2str(numDimsPerGroupCands), ...
  datestr(now,'ddmm-hhMMss') );

% Compute some statistics to help with the initialization
th = rand(1000, numDims); fth = func(th);
meanFth = mean(fth);
stdFth = std(fth);

% Parameters for additive Bayesian optimization
boParams.optPtStdThreshold = 0.002;
boParams.alBWLB = 1e-5;
boParams.alBWUB = 5;
boParams.numInitPts = 10; %20; % min(20, numDims);
boParams.commonNoise = 0.01 * stdFth;
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

% From here on customize each parameters separately.
% Known Decomposition
boKDParams = boParams;
boKDParams.decompStrategy = 'known';
boKDParams.diRectParams.maxevals = ceil(numDiRectEvals/trueNumGroups);
boKDParams.diRectParams.maxits = inf;
% Known grouping but not decomposition
boUDParams = boParams;
if numDims < 12
  boUDParams.decompStrategy = 'learn';
else
  boUDParams.decompStrategy = 'partialLearn';
end
boUDParams.diRectParams.maxevals = ceil(numDiRectEvals/trueNumGroups);
boUDParams.diRectParams.maxits = inf;
% The rest - arbitrary decompositions
boAddParams = boParams;
boAddParams.decompStrategy = 'partialLearn';
boAddParams.diRectParams.maxits = inf;

totalNumQueries = numIters + boParams.numInitPts;
% Initialize arrays for storing the history
boKDHistories = zeros(numExperiments, totalNumQueries);
boUDHistories = zeros(numExperiments, totalNumQueries);
boAddHistories = zeros(numExperiments, totalNumQueries, numdCands);
randHistories = zeros(numExperiments, totalNumQueries);
% For storing simple regret values
boKDSimpleRegrets = zeros(numExperiments, totalNumQueries);
boUDSimpleRegrets = zeros(numExperiments, totalNumQueries);
boAddSimpleRegrets = zeros(numExperiments, totalNumQueries, numdCands);
randSimpleRegrets = zeros(numExperiments, totalNumQueries);
% For storing cumulative regret values
boKDCumRegrets = zeros(numExperiments, totalNumQueries);
boUDCumRegrets = zeros(numExperiments, totalNumQueries);
boAddCumRegrets = zeros(numExperiments, totalNumQueries, numdCands);
randCumRegrets = zeros(numExperiments, totalNumQueries);

% First Call Direct
diRectOpts.maxevals = totalNumQueries;
diRectOpts.maxits = inf;
diRectOpts.showits = true;
[~, ~, diRectHist, ~, diRectHistory] = diRectWrap(func, bounds, diRectOpts);
[diRectSimpleRegret, diRectCumRegret] = getRegrets(trueMaxVal, diRectHistory);

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

%   % Learn Decomposition
%   fprintf('\nKnown Grouping Unknown Decomposition\n');
%   [decompUD, boUDParams] = ...
%     getDecompForParams(numDims, trueNumDimsPerGroup, boUDParams);
%   [~, ~, ~, valHistUD] = ...
%     bayesOptDecompAddGP(func, decompUD, bounds, numIters, boUDParams);
%   [sR, cR] = getRegrets(trueMaxVal, valHistUD);
%   boUDHistories(expIter, :) = valHistUD';
%   boUDSimpleRegrets(expIter, :) = sR';
%   boUDCumRegrets(expIter, :) = cR';

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
    [sR, cR] = getRegrets(trueMaxVal, valHistAdd);
    boAddHistories(expIter, :, candIter) = valHistAdd';
    boAddSimpleRegrets(expIter, :, candIter) = sR';
    boAddCumRegrets(expIter, :, candIter) = cR';
  end

  % Random
  randQueries = bsxfun(@plus, ...
    bsxfun(@times, rand(totalNumQueries, numDims), ...
      (bounds(:,2) - bounds(:,1))' ), bounds(:,1)' );
  randHistories(expIter, :) = func(randQueries)';
  [sR, cR] = getRegrets(trueMaxVal, randHistories(expIter, :)');
  randSimpleRegrets(expIter, :) = sR';
  randCumRegrets(expIter, :) = cR';

  % Save Results at each iteration
  save(saveFileName, 'numDims', 'trueNumDimsPerGroup', 'func', ...
    'funcProperties', 'trueMaxVal', 'numDimsPerGroupCands',  ...
    'numIters', 'totalNumQueries', 'numdCands', ...
    'boKDHistories', 'boUDHistories', 'boAddHistories', 'randHistories', ...
    'diRectHist', 'diRectHistory', ...
    'boKDSimpleRegrets', 'boUDSimpleRegrets', 'boAddSimpleRegrets', ...
    'randSimpleRegrets', 'diRectSimpleRegret', ...
    'boKDCumRegrets', 'boUDCumRegrets', 'boAddCumRegrets', 'randCumRegrets', ...
    'diRectCumRegret' );
  
end

  % Now Plot them out
  clearvars -except saveFileName;
  load(saveFileName);
  plotGPBResults;

