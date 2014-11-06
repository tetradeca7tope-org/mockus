function [maxVal, maxPt, boQueries, boVals, history] = bayesOptAddGP(...
  oracle, decomposition, bounds, numIters, params)
% Performs BO to maximize the function given in the oracle. Read bayesOptGP.m
% for explanation on the inputs and outputs. The only difference here is the new
% input argument decomposition. This is a matlab cell where each element is a
% vector containing the coordinates of the group.

  % Prelims
  numDims = size(bounds, 1);
  numGroups = numel(decomposition);
  dummyPts = zeros(0, numDims); % to build the GPs
  MAX_THRESHOLD_EXCEEDS = 5;
  NUM_ITERS_PER_PARAM_RELEARN = 20;

  % Do some diagnostics on the decomposition and print them out
  relevantCoords = cell2mat(decomposition');
  numRelevantCoords = numel(relevantCoords);
  if numRelevantCoords ~= numel(unique(relevantCoords))
    error('The Same coordinate cannot appear in different groups');
  end
  fprintf('# Groups: %d, %d/%d coordinates are relevant\n', ...
    numGroups, numRelevantCoords, numDims);

  % If Init Points are not given, initialize
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
    fprintf('WARNING: diRectParams not given. DiRect will use default vals.\n');
%     params.diRectParams.maxits = 8;
  end
  if ~isfield(params, 'useFixedBandWidth')
    params.useFixedBandWidth = false;
  end

  % Determine values to be sent to the GP
  % ============================================================================
  % The Common Mean Function
  % -----------------------------------------------------
    if isfield(params, 'commonMeanFunc')
      gpHyperParams.commonMeanFunc = params.commonMeanFunc;
    else
      if isfield(params, 'meanFuncValue')
        gpHyperParams.commonMeanFunc = ...
          @(arg) params.meanFuncValue * ones(size(arg,1), 1);
      else
        gpHyperParams.commonMeanFunc = [];
      end
    end
  % The Indiviual GP Mean Functions
  % -----------------------------------------------------
    if ~isfield(params, 'meanFuncs') | isempty(params.meanFuncs)
      gpHyperParams.meanFuncs = @(arg) zeros( size(arg, 1), 1);
    else
      gpHyperParams.meanFuncs = params.meanFuncs; 
      % if only one is given, they will be dealt to all GPs in
      % addGPMargLikelihood.
    end
  % The Common Noise Parameters
  % -----------------------------------------------------
    if ~isfield(params, 'commonNoise') | isempty(params.commonNoise)
      gpHyperParams.commonNoise = 0.01 * std(initVals);
    else
      gpHyperParams.commonNoise = params.commonNoise;
    end
  % The Indiviual Noise Functions
  % -----------------------------------------------------
    if ~isfield(params, 'noises') | isempty(params.noises)
      gpHyperParams.noises = 0;
    else
      gpHyperParams.noises = params.noises; 
    end
  % The Scale Parameters
  % -----------------------------------------------------
    if isfield(params, 'fixPr')
      gpHyperParams.fixPr = params.fixPr;
    else
      gpHyperParams.fixPr = false;
    end
    gpHyperParams.useSamePr = params.useSamePr;
    if params.useSamePr
      gpHyperParams.sigmaPrRange = params.sigmaPrRange;
    else
      gpHyperParams.sigmaPrRanges = params.sigmaPrRanges;
    end
  % The Bandwidth 
  % -----------------------------------------------------
    gpHyperParams.useSameSm = params.useSameSm;
    if params.useFixedBandWidth
      gpHyperParams.fixSm = true;
      if params.useSameSm, gpHyperParams.sigmaSmRange = params.sigmaSmRange;
      else, gpHyperParams.sigmaSmRanges = params.sigmaSmRanges;
      end
    else % This BO algorithm will set the bw via its own procedure
      alBWLB = params.alBWLB;
      alBWUB = params.alBWUB;
      % Set an initial bw. This will change as the algorithm progresses
      alCurrBW = alBWUB;
    end

  % Define the following before proceeding
  boQueries = [initPts; zeros(numIters, numDims)];
  boVals = [initVals; zeros(numIters, 1)];
  history = [max(initVals(cumsum(triu(ones(length(initVals))))))'; ...
             zeros(numIters, 1)];
  threshExceededCounter = 0;
  currMaxVal = -inf;
  
  fprintf('Peforming BO (dim = %d)\n', numDims);
  for boIter = 1:numIters

    if mod(boIter, 20) == 0
      fprintf('Additive GP BO iter %d/ %d\n', boIter, numIters);
    end

    % Prelims
    numBoPts = numInitPts + boIter - 1;
    currBoQueries = boQueries(1:numBoPts, :);
    currBoVals = boVals(1:numBoPts);

    % First redefine ranges for the GP bandwidth if needed
    if ~params.useFixedBandWidth & ...
       ( mod (boIter-1, NUM_ITERS_PER_PARAM_RELEARN) == 0 | ...
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
        gpHyperParams.sigmaSmRanges = alBWLB * ones(numGroups, 1);

      else
        gpHyperParams.fixSm = false;
        % Use same bandwidth for now. TODO: modify to allow different bandwidths
        gpHyperParams.useSameSm = true;
        gpHyperParams.sigmaSmRange = [alBWLB, alBWUB];
      end % alBWUB == alBWLB
      % Obtain the optimal GP parameters
      [~, ~, ~, ~, ~, ~, alCurrBWs, alCurrScales] = ...
        addGPMargLikelihood( currBoQueries, currBoVals, dummyPts, ...
          decomposition, gpHyperParams );
      alCurrBW = alCurrBWs(1); %TODO: modify to allow different bandwidths
      fprintf('Picked bw: %0.4f (%0.4f, %0.4f), Scale: %0.4f\n', ...
        alCurrBW, alBWLB, alBWUB, alCurrScales(1));

    end % if ~params.useFixedBandWidth ...

    % Now build the GP
    currGPParams.commonMeanFunc = gpHyperParams.commonMeanFunc;
    currGPParams.meanFuncs = gpHyperParams.meanFuncs;
    currGPParams.commonNoise = gpHyperParams.commonNoise;
    currGPParams.noises = gpHyperParams.noises;
    currGPParams.sigmaSms = alCurrBWs;
    currGPParams.sigmaPrs = alCurrScales;
    [~,~,~,~, combFuncH, funcHs] = addGPRegression(currBoQueries, ...
      currBoVals, dummyPts, decomposition, currGPParams);

    % Now obtain the next point
    [nextPt, ~, nextPtStd] = getNextQueryPt(params, combFuncH, funcHs, ...
      decomposition, currBoVals, bounds);
    % If it is too close, perturn it a bit
    if min( sqrt( sum( bsxfun(@minus, currBoQueries, nextPt).^2, 2) ))/...
          alCurrBW< 1e-10 
      while min( sqrt( sum( bsxfun(@minus, currBoQueries, nextPt).^2,2)))/...
                 alCurrBW < 1e-10
        nextPt = projectToRectangle( ...
          nextPt' + 0.1 * alCurrBW * randn(numDims, 1), bounds)';
      end
    end

    % Determine the current best point
    nextPtVal = oracle(nextPt);
    if nextPtVal > currMaxVal
      currMaxVal = nextPtVal;
      currMaxPt = nextPt;
    end
    boQueries(numInitPts + boIter, :) = nextPt;
    boVals(numInitPts + boIter) = nextPtVal;
%     fprintf('#: %d, maxVal: %0.5f, currVal: %0.5f\n', ...
%       boIter, currMaxVal, nextPtVal);

    % Check if nextPtStd is too small
    if nextPtStd < params.optPtStdThreshold
      threshExceededCounter = threshExceededCounter + 1;
    else
      threshExceededCounter = 0;
    end

    % Store the current best value
    history(boIter+numInitPts) = currMaxVal;

  end % end for boIter

  maxVal = currMaxVal;
  maxPt = currMaxPt;

end


function [nextPt, nextPtMean, nextPtStd, nextPtUtil] = ...
  getNextQueryPt(params, combFuncH, funcHs, decomposition, boVals, bounds)

  % So here's the plan here. We will define a utility function over each GP and
  % pick the point that maximizes this utility function
  
  numGroups = numel(decomposition);
  numDims = size(bounds, 1);

  nextPt = zeros(1, numDims);
  nextPtUtils = zeros(numGroups, 1);

  for k = 1:numGroups
    
    coords = decomposition{k};
    currBounds = bounds(coords, :);
    currGPFuncH = funcHs{k};

    if strcmp(params.utilityFunc, 'UCB')
      utility = @(t) getUCBUtility(t, currGPFuncH, size(boVals, 1)); 
    else
      utility = @(t) params.utilityFunc(t, currGPFuncH);
    end

    % Run DiRect
    [currNextPtUtil, currNextPt, hist] = ...
      diRectWrap(utility, currBounds, params.diRectParams);

    % store nextPt in relevant coordinates
    nextPt(coords) = currNextPt;
    nextPtUtils(k) = currNextPtUtil;

  end

  % Finally return the mean and the standard deviation at nextPt
  [nextPtMean, nextPtStd] = combFuncH(nextPt);

end

