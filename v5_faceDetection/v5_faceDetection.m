% Experiment Set up for Face Detection

close all;
clear all;
LIBPATH = '~/libs/kky-matlab/';
addpath([LIBPATH 'utils']); 
addpath([LIBPATH 'ancillary']); 
addpath([LIBPATH 'GPLibkky']); 
addpath ../BOLibkky/
addpath ../addGPLibkky/
addpath ../utils/
addpath VJ/
warning off;


% Problem parameters
numExperiments = 1;
numIters = 400;
numXTrain = 200;
numXTest = 10;
numDimsPerGroupCands = [22 1 2 5 11]';
doubleParams = false;

% Load the data
load VJ/vjdata
HaarCascade = GetHaarCasade('VJ/HaarCascades/haarcascade_frontalface_alt.mat');

% Other derived Parameters
if doubleParams, numDims = 44;
else numDims = 22; end
vjOptions.numTrain = numXTrain;
vjOptions.doubleParams = doubleParams;
numDiRectEvals = min(20000, max(1000*numDims, 500));
numdCands = numel(numDimsPerGroupCands);

% Get the function
vjBounds = [ 0 4; 0 20; 0 20; 5 25; 0 20; 10 30; 10 30; 10 30; ...
  10 30; 20 40; 30 40; 50 60; 40 60; 40 60; 50 70; 50 70; 55 75; ...
  70 80; 70 90; 80 100; 100 120; 100 120];
% vjBounds = [ 0 2; 0 10; 5 15; 10 20; 10 20; 15 25; 15 25; 20 30; ...
%   20 30; 30 40; 30 40; 45 55; 50 60; 45 55; 65 75; 60 70; 65 75; ...
%   75 85; 80 90; 85 95; 100 110; 100 110];
nominalParams = [0.8227 6.9566 9.4985 18.4130 15.3241 21.0106 23.9188 24.5279 ... 
  27.1534 34.5541 39.1073 50.6105 54.6201 50.1697 66.6691 67.6989 69.2288 ...
  79.2491 87.6960 90.2533 104.7492 105.7611];
normalizedNominalParams = getNormParams(nominalParams, vjBounds);
func = @(t) vjWrap(XTrain, YTrain, HaarCascade, t, vjBounds, vjOptions);
bounds = [zeros(numDims, 1) ones(numDims, 1)];

% Ancillary stuff
resultsDir = 'results/';
saveFileName = sprintf('%svj%d-%d-%s-%s.mat', resultsDir, numDims, numXTrain,...
  mat2str(numDimsPerGroupCands), datestr(now,'ddmm-hhMMss') );

% Parameters for additive Bayesian optimization
boParams.optPtStdThreshold = 0.002;
boParams.alBWLB = 1e-5;
boParams.alBWUB = 5;
boParams.numInitPts = 0; %20; % min(20, numDims);
boParams.commonNoise = 0.3;
boParams.utilityFunc = 'UCB';
boParams.meanFuncs = [];
boParams.commonMeanFunc = @(arg) zeros(size(arg, 1), 1);
% boParams.commonMeanFunc = []; %@(arg) meanFth * ones(size(arg, 1), 1);
boParams.useSamePr = true;
boParams.useSameSm = true;
boParams.fixPr = false;
boParams.fixSm = false;
boParams.sigmaPrRange = [5 25];
boParams.useFixedBandwidth = false;

% The rest - arbitrary decompositions
boAddParams = boParams;
boAddParams.decompStrategy = 'partialLearn';
boAddParams.diRectParams.maxits = inf;

totalNumQueries = numIters + boParams.numInitPts;
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
nominalVal = func(normalizedNominalParams);

% First Call Direct
fprintf('First Running DiRect\n============================================\n');
diRectOpts.maxevals = totalNumQueries;
diRectOpts.maxits = inf;
diRectOpts.showits = true;
tic;
[~, diRectOptPt, diRectHist] = diRectWrap(func, bounds, diRectOpts);
diRectTime = toc;
[diRectHistory, diRectMaxVals, diRectCumRewards] = ...
  getDiRectResults(diRectHist, -inf, totalNumQueries);
fprintf('DiRectOpt = %.4f\n', diRectMaxVals(end));

for expIter = 1:numExperiments

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
    boAddParamsCurr.diRectParams.maxevals = ceil(numDiRectEvals/numCurrGroups);
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

  % Random
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

  % Save Results at each iteration
  save(saveFileName, 'numDims', 'numDimsPerGroupCands', ...
    'numIters', 'totalNumQueries', 'numdCands', ...
    'boAddHistories', 'randHistories', 'diRectHistory', ...
    'boAddMaxVals', 'randMaxVals', 'diRectMaxVals', ...
    'boAddCumRewards', 'randCumRewards', 'diRectCumRewards' , ...
    'boAddQueryPts', 'randQueryPts', 'diRectOptPt', ...
    'boAddTimes', 'randTimes', 'diRectTime');
  
end

  % Now Plot them out
  clearvars -except saveFileName;
  load(saveFileName);
  plotVJResults;

