function [mu, KPost, Mus, KPosts, combinedFuncH, funcHs, ...
  sigmaSmOpts, sigmaPrOpts] = ...
  addGPUGMargLikelihood(X, y, Xtest, decomposition, hyperParams)
% The same as addGPMargLikelihood but here we don't know the grouping.

  % prelims
  numDims = size(X, 1);
  numPts = size(X, 2);
  useSameSm = hyperParams.useSameSm;
  useSamePr = hyperParams.useSamePr;
  diRectOptions.maxits = 8;

  if iscell(decomposition)
    findDecomposition = false;
    numOuterIters = 1;
  else
    findDecomposition = true;
    numAOptIters = hyperParams.numAOptIters;
    numAOptInits = hyperParams.numAOptInits;
    numOuterIters = 2;

    % Initialize the decomposition
    temp = decomposition;
    M = numel(temp);
    cnt = 0;
    decomposition = cell(M, 1);
    for i = 1:M
      decomposition{i} = (cnt+1):(cnt + temp(i));
      cnt = cnt + temp(i);
    end
    p = numel([decomposition{:}]);

    % Potential Init Points for A
    if isempty(hyperParams.potInitAs), potInitAs = zeros(numDims, p, 0);
    else, potInitAs = hyperParams.potInitAs;
    end
  end
  numGroups = numel(decomposition);

  % Set the hyperparameters for each GP
  % Common Mean Function 
  % -----------------------------------------
  if isempty(hyperParams.commonMeanFunc)
    commonMeanFunc = @(arg) mean(y) * ones(size(arg,1), 1);
  else
    commonMeanFunc = hyperParams.commonMeanFunc;
  end
  % Common Noise parameter
  % -------------------------------------
    commonNoise = hyperParams.commonNoise;
  % Mean Functions for each GP
  % -------------------------------------
    if isempty(hyperParams.meanFuncs)
      hyperParams.meanFuncs = @(arg) zeros( size(arg,1), 1);
    end
    if numel(hyperParams.meanFuncs) == 1
      meanFuncs = cell(numGroups, 1);
      [meanFuncs{:}] = deal(hyperParams.meanFuncs);
    else
      meanFuncs = hyperParams.meanFuncs;
    end
  % Noise Parameters
  % -------------------------------------
    if numel(hyperParams.noises) == 1
      noises = hyperParams.noises * ones(numGroups, 1);
    else
      noises = hyperParams.noises;
    end

  % Prepare bounds for the parameters. 
  % -------------------------------------
  % If using the same bw across all groups
  if useSameSm
    if isempty(hyperParams.sigmaSmRange)
      sigmaSmRange = norm(std(X));
    else
      sigmaSmRange = hyperParams.sigmaSmRange;
    end
    % Set the Bounds for Optimization
    if size(sigmaSmRange, 2) == 2, sigmaSmBound = sigmaSmRange;
    else, sigmaSmBound = sigmaSmRange * [0.01 100];
    end
  else
  % Different Bandwidths for the different groups
    if isempty(hyperParams.sigmaSmRanges)
      sigmaSmRanges = zeros(numGroups, 1);
      for k = 1:numGroups
        sigmaSmRanges(k) = norm( std( X(:, decomposition{k}) ) );
      end
    else
      sigmaSmRanges = hypeParams.sigmaSmRanges;
    end
    % Set the bounds for optimization
    if size(sigmaSmRanges, 2) == 2, sigmaSmBounds = sigmaSmRanges;
    else sigmaSmBounds = sigmaSmRanges * [0.01 100];
    end
  end

  % Now do the same for the scale parameter
  if useSamePr
    if isempty(hyperParams.sigmaPrRange)
      sigmaPrRange = std(y);
    else
      sigmaPrRange = hyperParams.sigmaPrRange;
    end
    % Set the Bounds for Optimization
    if size(sigmaPrRange, 2) == 2, sigmaPrBound = sigmaPrRange;
    else, sigmaPrBound = sigmaPrRange * [0.1 10];
    end
  else
  % Different Scales for the different groups
    if isempty(hyperParams.sigmaPrRanges)
      sigmaPrRanges = std(y) * ones(numGroups, 1);
    else
      sigmaPrRanges = hyperParams.sigmaPrRanges;
    end
    % Set the Bounds for Optimization
    if size(sigmaPrRanges, 2) == 2, sigmaPrBounds = sigmaPrRanges;
    else, sigmaPrBounds = sigmaPrRanges * [0.1 10];
    end
  end

  % Not adding noise here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Since there are several noise parameters, its best if the caller explicitly
  % sets the noise parametrs.

  % If fixing the smoothness/ scale, set the optimal values right away
  if hyperParams.fixSm
    if useSameSm, sigmaSmOpts = sigmaSmRange * ones(numGroups, 1);
    else, sigmaSmOpts = sigmaSmRange
    end
  end
  if hyperParams.fixPr
    if useSamePr, sigmaPrOpts = sigmaPrRange * ones(numGroups, 1);
    else, sigmaPrOpts = sigmaPrRange;
    end
  end

  % Prelims for the Optimization
  oneVec = ones(numGroups, 1);

  % ============================================================================
  % ============================================================================
  for outerIter = 1:numOuterIters

    % First optimize for the Smoothness and the Scale. We need to do this over
    % different scenarios.
    % ==========================================================================
    % 1. If fixing both parameters, nothing to be done
    % ------------------------------------------------
    if hyperParams.fixSm && hyperParams.fixPr

    % 2. If fixing only the Smoothness
    % ------------------------------------------------
    elseif hyperParams.fixSm
      % only optimize over the scale paraemeter 
      if useSamePr
        nlmlF = @(t) normMargLikelihood(sigmaSmOpts, t * oneVec, ...
          decomposition, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
        [~, optPr] = diRectWrap(nlmlF, sigmaPrBound, diRectOptions);
        sigmaPrOpts = optPr * oneVec;
      else
        nlmlF = @(t) normMargLikelihood(sigmaSmOpts, t, ...
          decomposition, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
        [~, sigmaPrOpts] = diRectWrap(nlmlF, sigmaPrBounds, diRectOptions);
      end

    % 3. If fixing only the Scale 
    % ------------------------------------------------
    elseif hyperParams.fixPr
      % only optimize over the bandwidth parameter
      if useSamePr
        nlmlF = @(t) normMargLikelihood(t * oneVec, sigmaPrOpts,  ...
          decomposition, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
        [~, optSm] = diRectWrap(nlmlF, sigmaSmBound, diRectOptions);
        sigmaSmOpts = optSm * oneVec;
      else
        nlmlF = @(t) normMargLikelihood(t, sigmaPrOpts, ...
          decomposition, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
        [~, sigmaSmOpts] = diRectWrap(nlmlF, sigmaSmBounds, diRectOptions);

      end

    % 4. If fixing neither
    % ------------------------------------------------
    else
      % 4.1 If using the same for both scale and bandwidth for all groups
      if useSameSm && useSamePr
        nlmlF = @(t) normMargLikelihood( t(1)*oneVec, t(2)*oneVec, ...
          decomposition, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
        diRectBounds = [sigmaSmBound; sigmaPrBound];
        [~, optParams] = diRectWrap(nlmlF, diRectBounds, diRectOptions);
        sigmaSmOpts = optParams(1) * oneVec;
        sigmaPrOpts = optParams(2) * oneVec;

      % 4.2 If using the same smoothness for all groups
      elseif useSameSm
        nlmlF = @(t) normMargLikelihood( t(1)*oneVec, t(2:end) * oneVec, ...
          decomposition, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
        diRectBounds = [sigmaSmBound; sigmaPrBounds];
        [~, optParams] = diRectWrap(nlmlF, diRectBounds, diRectOptions);
        sigmaSmOpts = optParams(1) * oneVec;
        sigmaPrOpts = optParams(2:end);

      % 4.3 If using the same scale
      elseif useSamePr
        nlmlF = @(t) normMargLikelihood( t(1:numGroups), t(end) * oneVec, ...
          decomposition, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
        diRectBounds = [sigmaSmBounds; sigmaPrBound];
        [~, optParams] = diRectWrap(nlmlF, diRectBounds, diRectOptions);
        sigmaSmOpts = optParams(1:numGroups); 
        sigmaPrOpts = optParams(end) * oneVec;

      % 4.4 If using different scales and bws
      else
        nlmlF = @(t) normMargLikelihood( t(1:numGroups), t(numGroups+1:end), ...
          decomposition, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
        diRectBounds = [sigmaSmBounds; sigmaPrBounds];
        [~, optParams] = diRectWrap(nlmlF, diRectBounds, diRectOptions);
        sigmaSmOpts = optParams(1:numGroups);
        sigmaPrOpts = optParams(numGroups+1 : end); 

      end %  if useSameSm && useSamePr

    end % if fixSm & fixPr
    % ==========================================================================
    % phew! that was crazy ...

    % Now optimize over the decomposition
    if findDecomposition & (outerIter ~= numOuterIters)
      negMLF = @(T) negNormRotMargLikelihood( sigmaSmOpts, sigmaPrOpts, ...
        decomposition, T, X, Y, meanFuncs, commonMeanFunc, noises, commonNoise);
      % TODO: Optimize w.r.t A here and then project this to the set of
      % permutation matrices.
      A = ce3AOptimize(negMLF, numDims, p, [], []);

    end

  end % for outerIter = 1:numOuterIters
  % ============================================================================
  % ============================================================================

  % Finally train the GP
  gpHyperParams.commonMeanFunc = commonMeanFunc;
  gpHyperParams.meanFuncs = meanFuncs;
  gpHyperParams.commonNoise = commonNoise;
  gpHyperParams.noises = noises;
  gpHyperParams.sigmaSms = sigmaSmOpts;
  gpHyperParams.sigmaPrs = sigmaPrOpts;
  [mu, KPost, Mus, KPosts, combinedFuncH, funcHs] = ...
    addGPRegression(X, y, Xtest, decomposition, gpHyperParams);

end


function nlml = normMargLikelihood(sigmaSms, sigmaPrs, decomposition, ...
  X, y, meanFuncs, commonMeanFunc, noises, commonNoise)
% Computes the normalized log likelihood of the model.

  % prelims
  numData = size(X, 1);
  Ky = combinedKernelNoise(X, decomposition, sigmaSms, sigmaPrs, noises, ...
                           commonNoise);
  L = stableCholesky(Ky);
  y_ = y - combinedMeanFunc(X, commonMeanFunc, meanFuncs, decomposition);
  alpha = L' \ (L \ y_);
  nlml = -1/2 * y_' * alpha - sum(log(diag(L))) - numData/2 * log(2*pi);
end

