function [maxVal, maxPt, boQueries, boVals, history] = bayesOptAddUDGP(...
  oracle, d, M,  bounds, numIters, params)
% Performs BO to maximize the function given in the oracle. This is like
% bayesOptAddGP but now we do not know the decomposition. So at each iteration
% we will attempt to maximize the marginal likelihood to find the decomposition.
% The decomposition is treated as a parameter of the kernel and is maximized
% along with the bandwidth and the scale.
% Inputs
% d: the dimensionality of an additive group.
% M: the number of such groups. Otherwise, a vector with the dimensionality for
%    each group.

  % Prelims
  numDims = size(bounds, 1);
  D = numDims;
  MAX_THRESHOLD_EXCEEDS = 5;
  NUM_ITERS_PER_PARAM_RELEARN = 20;
  
  if isscalar(M)
    decomposition = cell(M,1);
    for i = 1:M, decomposition{i} = ( (i-1)*d + 1): (i*d); end
  else
    dimsPerGroup = M;
    M = numel(dimsPerGroup);
    decomposition = cell(M, 1);
    cnt = 0;
    for i = 1:M
      decomposition{i} = (cnt+1):(cnt+dimsPerGroup(i));
      cnt = cnt + dimsPerGroup(i);
    end
  end
  numGroups = numel(decomposition); % The same as M
  p = numel([decomposition{:}]);

  % Print out some diagnostics
  relevantCoords = cell2mat(decomposition');
  numRelevantCoords = numel(relevantCoords);
  fprintf('# Groups: %d, %d/%d coordinates are relevant\n', numGroups, ...
    numRelevantCoords, numDims);

  % If Init points are not given, initialize
  if ~isfield(params, 'initPts') | isempty(params.initPts)
    initPts = boGetInitPts(bounds, params.numInitPts);
    initVals = oracle(initPts);
  else
    initPts = params.initPts;
    initVals = params.initVals;
  end
  numInitPts = size(initPts, 1);

  % Check for parameters expected in params
  if ~isfield(params, 'diRectParams')
    params.diRectParams = struct;
    fprintf('WARNING: diRectParams not given. DiRect will use default vals\n');
  end
  if ~isfield(params, 'useFixedBandWidth')
    params.useFixedBandWidth = false;
  end

  % Determine values to be sent to the GP
  % ============================================================================
  % The Common Mean Function
  % --------------------------------------------------------------
    if isfield(params, 'commonMeanFunc')
      gpHyperParams.commonMeanFunc = params.commonMeanFunc;
      % NOTE: the commonMeanFunc's arguments are the variables in the original
      % space (before the transformation).
    else
      if isfield(params, 'meanFuncValue')
        gpHyperParams.commonMeanFunc = ...
          @(arg) params.meanFuncValue * ones( size(arg,1), 1);
      else
        gpHyperParams.commonMeanFunc = [];
      end
    end
  % The Individual GP Mean Functions
  % --------------------------------------------------------------
    if ~isfield(params, 'meanFuncs') | isempty(params.meanFuncs)
      gpHyperParams.meanFuncs = @(arg) zeros( size(arg, 1), 1);
    else
      gpHyperParams.meanFuncs = params.meanFuncs;
    end
  % The Common noise parameters
  % --------------------------------------------------------------
    if ~isfield(params, 'commonNoise') | isempty(params.commonNoise)
      gpHyperParams.commonNoise = 0.01 * std(initVals);
    else
      gpHyperParams.commonNoise = params.commonNoise;
    end
  % The Individual noise params
  % --------------------------------------------------------------
    if ~isfield(params, 'noises') | isempty(params.noises)
      gpHyperParams.noises = 0;
    else
      gpHyperParams.noises = params.noises;
    end
  % Currently my GP implementation uses the same bw and scale for all kernels
  % so need for parameters such as useSamePr, useSameSm
  if isfield(params, 'fixPr')
    gpHyperParams.fixPr = params.fixPr;
  else
    gpHyperParams.fixPr = false;
  end
  gpHyperParams.fixSm = params.useFixedBandWidth;
  if params.fixSm
    % Then we should fix the bandwidth and continue.
    gpHyperParams.sigmaSmRange = params.sigmaSmRange;
  else
    alBWLB = params.alBWLB;
    alBWUB = params.alBWUB;
    alCurrBW = alBWUB;
  end

  % Define the following before proceeding
  dummyPts = zeros(0, numDims);
  boQueries = [initPts; zeros(numIters, numDims)];
  boVals = [initVals; zeros(numIters, 1)];
  history = [max(initVals(cumsum(triu(ones(length(initVals))))))'; ...
             zeros(numIters, 1)];
  threshExceededCounter = 0;
  currMaxVal = max(boVals(1:numInitPts));
  % Other GP Hyper params
  gpHyperParams.sigmaPrRange = [];

  fprintf('Peforming BO (dim = %d)\n', numDims);
  for boIter = 1:numIters

    if mod(boIter, 20) == 0
      fprintf('Additive GP BO iter %d/ %d\n', boIter, numIters);
    end

    % Prelims
    numBoPts = numInitPts + boIter - 1;
    currBoQueries = boQueries(1:numBoPts, :);
    currBoVals = boVals(1:numBoPts);

    % First redefine ranges for the GP bw if needed
    if ~params.useFixedBandWidth & ...
      ( mod(boIter-1, NUM_ITERS_PER_PARAM_RELEARN) == 0 | ...
        threshExceededCounter == MAX_THRESHOLD_EXCEEDS )

      if threshExceededCounter == MAX_THRESHOLD_EXCEEDS
        alBWUB = max(alBWLB, 0.9 * alCurrBW);
        threshExceededCounter = 0;
        fprintf('Threshold Exceeded %d times - Reducing BW\n', ...
          MAX_THRESHOLD_EXCEEDS);
      else
        alBWUB = max(alBWLB, 0.9 * alBWUB);
      end

      % Define the BW range for addGPMargLikelihood
      if alBWUB == alBWLB
        gpHyperParams.fixSm = true;
        gpHyperParams.useSameSm = true;

      else
        gpHyperParams.fixSm = true;
        gpHyperParams.sigmaSmRange = [alBWLB, alBWUB];
      end

      % Now obtain the optimal hyper parameters
      [~, ~, ~, ~, ~, ~, ~, alCurrBWs, alCurrScales, A, decomposition] = ...
        addGPRotMargLikelihood(currBoQueries, currBoVals, dummyPts, d, M, ...
                               gpHyperParams);
      alCurrBW = alCurrBWs(1);
      fprintf('Picked BW: %0.4f (%0.4f, %0.4f), Scale: %0.4f\n', ...
        alCurrBW, alBWLB, alBWUB, alCurrScales(1));
      A,

    end % if ~params.useFixedBandWidth

    % Now build the GP
    currGPParams.commonMeanFunc = gpHyperParams.commonMeanFunc;
    currGPParams.meanFuncs = gpHyperParams.meanFuncs;
    currGPParams.commonNoise = gpHyperParams.commonNoise;
    currGPParams.noises = gpHyperParams.noises;
    currGPParams.sigmaSms = alCurrBWs;
    currGPParams.sigmaPrs = alCurrScales;
    Z = currBoQueries * A;
    [~,~,~,~, combZFuncH, funcHs] = addGPRegression(Z, ...
      currBoVals, dummyPts, decomposition, currGPParams);

    % Now obtain the next point
    [nextPt, ~, nextPtStd] = getNextQueryPt(A, params, combZFuncH, funcHs, ...
      decomposition, currBoVals, bounds);
    % If it is too close to any of the existing points, perturb it a bit
    if min( sqrt( sum( bsxfun(@minus, currBoQueries, nextPt).^2, 2) ))/...
          alCurrBW< 1e-10
      while min( sqrt( sum( bsxfun(@minus, currBoQueries, nextPt).^2,2)))/...
                 alCurrBW < 1e-10
        nextPt = projectToRectangle( ...
          nextPt' + 0.1 * alCurrBW * randn(numDims, 1), bounds)';
      end
    end

    % Evaluate and determine the current best point
    nextPtVal = oracle(nextPt);
    if nextPtVal > currMaxVal
      currMaxVal = nextPtVal;
      currMaxPt = nextPt;
    end
    boQueries(numInitPts + boIter, :) = nextPt;
    boVals(numInitPts + boIter) = nextPtVal;

    % Check if nextPtStd is too small
    if nextPtStd < params.optPtStdThreshold
      threshExceededCounter = threshExceededCounter + 1;
    else
      threshExceededCounter = 0;
    end

    % Store the current best value
    history(boIter+numInitPts) = currMaxVal;

  end % for boIter

  maxVal = currMaxVal;
  maxPt = currMaxPt;
  [boQueries, boVals],
  

end


function [nextPt, nextPtMean, nextPtStd, nextPtUtil] = ...
  getNextQueryPt(A, params, combZFuncH, funcHs, decomposition, boVals, bounds)

  numGroups = numel(decomposition);
  numDims = size(bounds, 1);

  nextZPt = zeros(1, numDims);
  nextPtUtils = zeros(numGroups, 1);

  groupOrder = randperm(numGroups);
  for k = groupOrder
    coords = decomposition{k};
    currBounds = bounds(coords, :);
    numGroupDims = size(currBounds, 1);
    currGPFuncH = funcHs{k};

    if strcmp(params.utilityFunc, 'UCB')
      utility = @(t) getUCBUtility(t, currGPFuncH, size(boVals, 1)); 
    else
      utility = @(t) params.utilityFunc(t, currGPFuncH);
    end

    % Optimize over this group by searching exhaustively.
    numCandidates = min(30^numGroupDims, 20000);
    groupCandidates = getUniformPtsInBounds(numCandidates, bounds);
    % Now transform them to the Z space of this group
    groupZCands = groupCandidates * A(:, coords)';

    % Run DiRect
    [currNextPtUtil, currNextPt, hist] = ...
      diRectWrap(utility, currBounds, params.diRectParams);
%     fprintf('%d, currNextPt: %s\n', k, mat2str(currNextPt));

    % store nextPt in relevant coordinates
    nextZPt(coords) = currNextPt;
    nextPtUtils(k) = currNextPtUtil;

  end

  % Finally get it in the X coordinates
  nextPt = nextZPt * A';

  [nextPtMean, nextPtStd] = combZFuncH(nextZPt);

end

