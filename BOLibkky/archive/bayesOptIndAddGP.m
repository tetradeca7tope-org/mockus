function [maxVal, maxPt, boQueries, boVals, history] = bayesOptAddGP(...
  oracle, decomposition, bounds, numIters, params)
% Performs BO on each point separately.

  % Prelims
  numDims = size(bounds, 1);
  numGroups = numel(decomposition);
  dummyPts = zeros(0, numDims);
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

  % Determine the Random Point in R^d first
  anchorPt = 0.5 * rand(numDims, 1);
  % For saving the best init point
  bestInitVals = zeros(numGroups, 1);
  bestInitPts = cell(numGroups, 1);

  if ~isfield(params, 'initPts') | isempty(params.initPts)
  % If Init Points are not given, initialize
    initPts = cell(numGroups, 1); completeInitPts = cell(numGroups, 1);
    initVals = cell(numGroups, 1);
    for k = 1:numGroups
      coords = decomposition{k};
      groupInitPts = boGetInitPts(bounds(coords, :), params.numInitPtsPerGroup);
      initPts{k} = groupInitPts;
      [vals, currInitQueryPts] = ...
        queryAtThisGroup(oracle, anchorPt, coords, groupInitPts);
      initVals{k} = vals;
      completeInitPts{k} = currInitQueryPts;
      [bestVal, bestIdx] = max(vals);
      bestInitVals(k) = bestVal;
      bestInitPts{k} = initPts{k}(bestIdx, :);
%       initPts{k}, initVals{k}, currInitQueryPts,
    end
  else
    initPts = params.initPts;
    initVals = params.initVals;
  end
  allInitPts = cell2mat(completeInitPts);
  allInitVals = cell2mat(initVals);
  numAllInitPts = size(initPts, 1);
  [bestInitVal, bestInitIdx] = max(allInitVals);

  % Check for parameters expected in params
  if ~isfield(params, 'diRectParams')
    params.diRectParams = struct;
    fprintf('WARNING: diRectParams not given. DiRect will use default vals.\n');
%     params.diRectParams.maxits = 8;
  end

  % Determine the values to be sent to the GP
  % ============================================================================
  % The Common Mean Function
  % -------------------------------------------------------------
%     if ~isfield(params, 'commonMeanFunc')
%       gpHyperParams.commonMeanFunc = commonMeanFunc;
%     else
%       if isfield(params, 'meanFuncValue')
%         gpHyperParams.commonMeanFunc = ...
%           @(arg) params.commonMeanFunc * ones(size(arg,1), 1);
%       else
%         gpHyperParams.commonMeanFunc = ...
%           @(arg) mean(Y) * ones(size(arg,1), 1);
%       end
%     end
  % Individual GP Mean Functions
  % -------------------------------------------------------------
    % This kind of gets hairy.

  % The Common Noise parameter
  % -------------------------------------------------------------
    if ~isfield(params, 'noise') | isempty(params.noise)
      gpHyperParams.noise = 0.01 * std(initVals{1});
    else
      gpHyperParams.noise = params.noise;
    end
  % The Individual Noise Functions
  % -------------------------------------------------------------
    if ~isfield(params, 'noises') | isempty(params.noises)
      gpHyperParams.noises = 0;
    else
      gpHyperParams.noises = params.noises;
    end
  % A Scale Range
  % -------------------------------------------------------------
    if isempty(params.sigmaPrRange)
      gpHyperParams.sigmaPrRange = std(allInitVals);
    else
      gpHyperParams.sigmaPrRange = params.sigmaPrRange;
    end

  % The bandwidth parameters
  alBWLBs = params.alBWLB * ones(numGroups, 1);
  alBWUBs = params.alBWUB * ones(numGroups, 1);
  alCurrBWs = alBWUBs;

  % Define the following before proceeding
  boQueries = [allInitPts; zeros(numIters, numDims)];
  boVals = [allInitVals; zeros(numIters, 1)];
  boGroupQueries = initPts;
  boGroupVals = initVals;
  currCombBestPt = zeros(numDims, 1);
  groupBestPts = cell(numGroups, 1);
  groupBestVals = bestInitVals;
  threshExceededCounter = zeros(numGroups, 1);
  bestPtAtEachIter = zeros(numIters, numDims);
  numGroupPts = params.numInitPtsPerGroup * ones(numGroups, 1);

  fprintf('Performing BO (dim = %d)\n', numDims);
  for boIter = 1:numIters

    if mod(boIter, 20) == 0
      fprintf('Ind-Add GP BO iter %d/ %d\n', boIter, numIters);
    end
    % Determine the current group
    currGroup = mod(boIter, numGroups);
    if currGroup == 0, currGroup = numGroups; end
    numGroupDims = numel(decomposition{currGroup});
    coords = decomposition{currGroup};

    currGpPts = numGroupPts(currGroup);

    if (mod(currGpPts-params.numInitPtsPerGroup, ...
            NUM_ITERS_PER_PARAM_RELEARN) == 0 | ...
      threshExceededCounter(currGroup) == MAX_THRESHOLD_EXCEEDS ) 
    % Then relearn the paramters

      if threshExceededCounter(currGroup) == MAX_THRESHOLD_EXCEEDS
        alBWUBs(currGroup) = max(1.05*alBWLB(currGroup), ...
                                  0.9*alCurrBWs(currGroup));
        threshExceededCounter(currGroup) = 0;
        fprintf('Threshold Exceeded %d times - Reducing BW\n', ...
                 MAX_THRESHOLD_EXCEEDS);
      else
        alBWUBs(currGroup) = max(1.05 * alBWLBs(currGroup), ...
                                  0.9 * alBWUBs(currGroup));
      end

      % Define BW Range for GP Marginal Likelihood
      dummyPts = zeros(0, numel(decomposition{currGroup}));
      gpHyperParams.sigmaPr = 0;
      gpHyperParams.sigmaSm = 0;
      gpHyperParams.meanFunc = params.meanFunc;
      gpHyperParams.sigmaSmRange = [alBWLBs(currGroup) alBWUBs(currGroup)];
      [~,~,funcH, optBw, optScale] = GPMargLikelihood( ...
        boGroupQueries{currGroup}, boGroupVals{currGroup}, dummyPts, ...
        gpHyperParams);
%       fprintf('Chosen BW/Scale: %0.4f, %0.4f\n', optBw, optScale);
      alCurrBWs(currGroup) = optBw;
      alCurrScales(currGroup) = optScale;

    end

    % Now Build the GP
    runTimeParams.retFunc = true;
    currGPParams.meanFunc = gpHyperParams.meanFunc;
    currGPParams.sigmaSm = alCurrBWs(currGroup);
    currGPParams.sigmaPr = alCurrScales(currGroup);
    currGPParams.noise = gpHyperParams.noise;
    dummyPts = zeros(0, numel(decomposition{currGroup}));
%     boGroupQueries{currGroup}, boGroupVals{currGroup}, dummyPts,
    [~, ~, ~, gpFuncH] = GPRegression(boGroupQueries{currGroup}, ...
      boGroupVals{currGroup}, dummyPts, currGPParams, runTimeParams);

    % Obtain the next Query Point
    [nextPt, ~, nextPtStd] = getNextQueryPt(params, gpFuncH, ...
      boGroupVals{currGroup}, bounds(coords, :));
    if min( sqrt( sum( bsxfun(@minus, boGroupQueries{currGroup}, ...
                              nextPt).^2, 2) ))/alCurrBWs(currGroup)< 1e-10 
      while min( sqrt( sum( bsxfun(@minus, boGroupQueries{currGroup}, ...
                                nextPt).^2, 2)))/alCurrBWs(currGroup)< 1e-10 
        nextPt = projectToRectangle( ...
          nextPt' + 0.1 * alCurrBWs(currGroup) * rand(numGroupDims, 1), ...
                bounds(coords, :))';
      end
    end

    [nextPtVal, nextQueryPt] = queryAtThisGroup(oracle, anchorPt, coords, ...
                                nextPt);
    boQueries(numAllInitPts+boIter, :) = nextQueryPt;
    boVals(numAllInitPts+boIter, :) = nextPtVal;
    boGroupQueries{currGroup} = [boGroupQueries{currGroup}; nextPt];
    boGroupVals{currGroup} = [boGroupVals{currGroup}; nextPtVal];

    % Store the best point in the group.
    if nextPtVal > groupBestVals(currGroup)
      groupBestVals(currGroup) = nextPtVal;
      groupBestPts{currGroup} = nextPt';
      currCombBestPt(coords,:) = nextPt';
    end
    bestPtAtEachIter(boIter,:) = currCombBestPt';

    if nextPtStd < params.optPtStdThreshold
      threshExceededCounter(currGroup)  = threshExceededCounter(currGroup) + 1;
    else
      threshExceededCounter(currGroup) = 0;
    end

    numGroupPts(currGroup) = numGroupPts(currGroup) + 1;

  end

  history = [max(allInitVals(cumsum(triu(ones(length(allInitVals))))))'; ...
             oracle(bestPtAtEachIter)];
  maxVal = history(end);
  maxPt = bestPtAtEachIter(end,:);
  fprintf('bayesOptIndAddGP val : %.4f\n', maxVal);

end


function [vals, queryPts] = queryAtThisGroup(oracle, anchorPt, coords, groupPts)

  numPts = size(groupPts, 1);
  queryPts = repmat(anchorPt', numPts, 1);
  queryPts(:, coords) = groupPts;
  vals = oracle(queryPts);

end


function [nextPt, nextPtMean, nextPtStd, nextPtUtil] = ...
  getNextQueryPt(params, gpFuncH, boVals, bounds)

  % First the utility function
  if strcmp(params.utilityFunc, 'EI')
    utility = @(t) getEIUtility(t, gpFuncH, max(boVals));
  elseif strcmp(params.utilityFunc, 'UCB')
    utility = @(t) getUCBUtility(t, gpFuncH, size(boVals, 1));
  else
    utility = @(t) params.utilityFunc(t, gpFuncH);
  end

  % Now run DiRect
  [nextPtUtil, nextPt, hist] = diRectWrap(utility, bounds, ...
    params.diRectParams);

  % Get Mean and Variance of Next Pt
  [nextPtMean, nextPtStd] = gpFuncH(nextPt);
end

