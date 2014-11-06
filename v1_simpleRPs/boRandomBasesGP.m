function [maxVal, maxPt, boQueries, boVals, history] = boRandomBaseGP(oracle, ...
  bounds, numIters, params)

  % Prelims
  numDims = size(bounds, 1);
  dummyPts = zeros(0, numDims); % to build GPs
  MAX_THRESHOLD_EXCEEDS = 5; %# times the variance is allowed to fall low
  NUM_ITERS_PER_PARAM_RELEARN = 20;

  % If Init points are not given, initialize
  if ~isfield(params, 'initPts') | isempty(params.initPts)
    if params.numInitPts == 1
      initPts = [ (bounds(:,2) + bounds(:,1))/2 ];
    else
      % Initialize k = params.numInitPts as follows: One point at the centre of
      % the grid. The other k-1 pts are chosen via k-means++ from a uniformly
      % sampled 1000 points
      if params.numInitPts < 10, numRandPts = 200;
      else numRandPts = 10 * params.numInitPts;
      end
      initPtCandidates = bsxfun(@plus, ...
        bsxfun(@times, rand(numRandPts, numDims), ...
               (bounds(:,2) - bounds(:,1))' ), bounds(:,1)');
      [~, initPts] = kmeanspp(initPtCandidates, params.numInitPts);
    end
    initVals = oracle(initPts);
  else
    initPts = params.initPts;
    initVals = params.initVals;
  end
  [initPts, initVals],
  numInitPts = size(initPts, 1);

  % Check for parameters expected in params
  if ~isfield(params, 'diRectParams')
    params.diRectParams = struct();
  end
  if ~isfield(params, 'useFixedBandwidth')
    params.useFixedBandwidth = false;
  end

  % Determine values to be sent to the GP
  % =========================================================
  % The Mean Function
  % --------------------------------------
    if ~isfield(params, 'meanFuncValue') | isempty(params.meanFuncValue)
    % If given, use it. Otherwise, GPLib will use the mean of the regressands
      gpHyperParams.meanFunc = [];
    else
      gpHyperParams.meanFunc = @(arg) params.meanFuncValue;
    end
  % The Scale
  % --------------------------------------
    if isfield(params, 'sigmaPr')
      % If the parameters specify the scale to use, then don't look beyond
      alCurrScale = params.sigmaPr;
    elseif isfield(params, 'scaleRange')
      alCurrScale = std(initVals);
      gpHyperParams.sigmaPrRange = params.scaleRange;
      gpHyperParams.sigmaPr =  0; % if we set this
    else ~isfield(params, 'scaleRange')
      % If the scale isn't given, but a range is given then we need to do
      % MargLikelihood. Specify this in the params.
      % If not given, then GPLib will use the regressands to estimate this.
      alCurrScale = std(initVals);
      gpHyperParams.sigmaPrRange = []; 
      gpHyperParams.sigmaPr =  0; % if we set this
    end
  % The Bandwidth
  % --------------------------------------
    if params.useFixedBandwidth
    % If use a fixed bandwidth, set it now. Otherwise we will set at each iter.
      gpHyperParams.sigmaSm = params.alBandWidth;
    else
      alBWLB = params.alBWLB;
      alBWUB = params.alBWUB;
      % Set an initial BW. This will change as the algorithm progresses.
      alCurrBW = alBWUB; 
    end
  % The Noise Level 
  % --------------------------------------
  if ~isfield(params, 'gpNoiseLevel')
    gpHyperParams.noise = 0.01 * std(initVals);
  else
    gpHyperParams.noise = params.gpNoiseLevel;
  end

  % Define the following before proceeding
  boQueries = [initPts; zeros(numIters, numDims)];
  boVals = [initVals; zeros(numIters, 1)];
  history = [max(initVals) * ones(numInitPts, 1); zeros(numIters, 1)];
  threshExceededCounter = 0;
  
  fprintf('Performing BO (dim = %d)\n', numDims);
  for boIter = 1:numIters

    if mod(boIter, 20) == 0
      fprintf('Bayesian Optimization iter %d/ %d\n', boIter, numIters);
    end
    % prelims
    numBoPts = numInitPts + boIter - 1;
    currBoQueries = boQueries(1:numBoPts, :);
    currBoVals = boVals(1:numBoPts);

    % First rebuild the GP if needed
    if ~params.useFixedBandwidth
      if mod(boIter-1, NUM_ITERS_PER_PARAM_RELEARN) == 0 | ...
          threshExceededCounter == MAX_THRESHOLD_EXCEEDS
        if threshExceededCounter == MAX_THRESHOLD_EXCEEDS
          alBWUB = max(alBWLB, 0.9 * alCurrBW);
          threshExceededCounter = 0;
          fprintf('Threshold exceeded %d times- Reducing Bandwidth\n', ...
                   MAX_THRESHOLD_EXCEEDS);
        end

        % Learn the optimal parameters for the GP Model.
        if alBWUB == alBWLB
          gpHyperParams.sigmaSm = alBWLB;
        else
          gpHyperParams.sigmaSm = 0;
          gpHyperParams.sigmaSmRange = [alBWLB, alBWUB];
        end
        [~, ~, ~, alCurrBW, alCurrScale] = GPMargLikelihood(currBoQueries, ...
          currBoVals, dummyPts, gpHyperParams);        
        fprintf('Picked bw, scale = %f , %f\n', alCurrBW, alCurrScale);
      end
    end

    % Build the GP: returns the function handle
    runTimeParams.retFunc = true;
    currGPParams.meanFunc = gpHyperParams.meanFunc;
    currGPParams.sigmaSm = alCurrBW;
    currGPParams.sigmaPr = alCurrScale;
    currGPParams.noise = gpHyperParams.noise;
    [~, ~, ~, gpFuncH] = GPRegression(currBoQueries, currBoVals, dummyPts, ...
                           currGPParams, runTimeParams);

    % Obtain the next query point and query
    % First prep the candidates
    NUM_AL_CANDIDATES = 100;
    A = rand(numDims, numDims);
    [A,~,~] = svd(A); % Just get the orthogonal directions
    dimPts = zeros(numDims, numDims);

    for dimIter = 1:numDims
    % Now, we will go through each direction and pick the points iteratively.
      a = A(:, dimIter);
      alCandidates = rand(NUM_AL_CANDIDATES, numDims);
      alCandidates = alCandidates * a * a' / (a'*a);
    
      [candMeans, candStds] = gpFuncH(alCandidates);
      candUtils = truncatedGaussianMean(candMeans, candStds, max(currBoVals));
      [~, nextPtIdx] = max(candUtils);
      nextPtInDim = alCandidates(nextPtIdx, :);
      dimPts(:, dimIter) = nextPtInDim;
    end
    nextPt = sum(dimPts, 2);
    nextPt = projectToRectangle(nextPt, bounds)';
    [~, nextPtStd] = gpFuncH(nextPt);

    % Query next pt 
    [nextPtVal] = oracle(nextPt);
    fprintf('#: %d, Value: %0.5f, Pt: %s\n', boIter,nextPtVal, mat2str(nextPt));
    boQueries(numInitPts + boIter, :) = nextPt;
    boVals(numInitPts + boIter) = nextPtVal;

    % Check if nextPtStd is too small
    if nextPtStd < params.optPtStdThreshold
      threshExceededCounter = threshExceededCounter + 1;
    else
      threshExceededCounter = 0;
    end

    % Store the current best value
    history(boIter+numInitPts) = max(history(boIter+numInitPts-1), nextPtVal);

  end % for boIter

  [maxVal, maxIdx] = max(boVals);
  maxPt = boQueries(maxIdx, :);
    
end


function [nextPt, nextPtMean, nextPtStd, nextPtUtil] = ...
  getNextPt(params, gpFuncH, boVals, bounds)
% This is what this function should do. It should maximize the utility 

  % First create the utility function to be maximized
  if strcmp(params.utilityFunc, 'EI')
    utility = @(t) getEIUtility(t, gpFuncH, max(boVals));
  else
    utility = @(t) params.utilityFunc(t, gpFuncH);
  end

  % Now run DiRect
  [nextPtUtil, nextPt, hist] = diRectWrap(utility, bounds, params.diRectParams);

  % Now get the mean and variance of nextPt
  [nextPtMean, nextPtStd] = gpFuncH(nextPt);

end

% Utility Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Expected Improvement
function util = getEIUtility(x, gpFuncH, trunc)
  [mu, s] = gpFuncH(x);
  util = truncatedGaussianMean(mu, s, trunc);
end

