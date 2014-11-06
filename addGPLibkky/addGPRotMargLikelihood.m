function [mu, KPost, Mus, KPosts, combinedXFuncH, combinedZFuncH, funcHs, ...
  sigmaSmOpts, sigmaPrOpts, A, learnedDecomp] = ...
  addGPRotMargLikelihood(X, y, Xtest, d, M, hyperParams)
% This function does the following. It attempts to find a suitable rotation
% matrix A \in D x p where p= dM so that Z = A'X forms a good basis for the
% additive function. We will then use an additive GP on the components of Z.
% Inputs
%   X, y, Xtest: respectively training inputs, outputs and testing data.
%   d: the dimensionality of an additive groups.
%   M: the number of such groups. otherwise a vector with the number in each
%      group.
%   hyperParams: see other GP implementations

  % prelims
  D = size(X, 2);
  numDims = D;
  numPts = size(X, 1);
  % for diRect
  diRectOptions.maxits = 8;

  if ~isfield(hyperParams, 'numOuterInits') | isempty(hyperParams.numOuterInits)
    numOuterInits = 2;
  else numOuterInits = hyperParams.numOuterInits;
  end
  if ~isfield(hyperParams, 'numBwSigmaDiRectIters') | ...
    isempty(hyperParams.numBwSigmaDiRectIters)
    numBwSigmaDiRectIters = 8;
  else numBwSigmaDiRectIters = hyperParams.numBwSigmaDiRectIters;
  end

  % Define the decomposition
  if isscalar(M)
    decomposition = cell(M, 1);
    for i = 1:M
      decomposition{i} = ((i-1)*d + 1): (i*d); end
  else
    dimsPerGroup = M;
    M = numel(dimsPerGroup);
    decomposition = cell(M, 1);
    cnt = 0;
    for i = 1:M
      decomposition{i} = (cnt+1):(cnt + dimsPerGroup(i));
      cnt = cnt + dimsPerGroup(i);
    end
  end
  p = numel([decomposition{:}]);
  numGroups = numel(decomposition);
     
  % Set the Hyperparameters for each GP
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

  % Initialization
  for initIter = 1:numOuterInits
    
    A = randn(D, p); A = orth(A);
    Z = X * A;

    % Define Bounds for optimization of the bandwidth and scale
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % bandwidth 
    if isempty(hyperParams.sigmaSmRange)
    % TODO: work with log bounds for sigma ??
      Zk = Z(:, decomposition{1}(:));
      sigmaSmBound = norm(std(Zk)) * [0.01 100];
    else
      sigmaSmRange = hyperParams.sigmaSmRange;
      if size(sigmaSmRange, 2) == 2, sigmaSmBound = sigmaSmRange;
      else sigmaSmBound = [0.01 100] * sigmaSmRange;
      end
    end
    % Scale
    if isempty(hyperParams.sigmaPrRange)
      sigmaPrBound = std(y) * [0.1 10];
    else
      sigmaPrRange = hyperParams.sigmaPrRange;
      if size(sigmaPrRange, 2) == 1, sigmaPrBound = sigmaPrRange;
      else sigmaPrBound = sigmaPrRange * [0.1 10];
      end
    end

    diRectOptions.maxits = numBwSigmaDiRectIters;

    oneVec = ones(numGroups, 1);
    nlmlF = @(t) normRotMargLikelihood( t(1)*oneVec, t(2)*oneVec, ...
      decomposition, A, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
    diRectBounds = [sigmaSmBound; sigmaPrBound];
    [~, optParams] = diRectWrap(nlmlF, diRectBounds, diRectOptions);
    sigmaSmOpts = optParams(1) * oneVec;
    sigmaPrOpts = optParams(2) * oneVec;

    % Now optimze w.r.t A
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    negMLF = @(T) negNormRotMargLikelihood( sigmaSmOpts, sigmaPrOpts, ...
      decomposition, T, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
%     A = ce3AOptimize(negMLF, D, p, true);
%     A = ce4AOptimize(negMLF, D, p, [], []);
    % If optimizing over just permutation matrices
    A = decompOptimize(negMLF, D, p, [], []);
    [~, permutedOrder] = orthToPermutation(A);
    
  end % for initIter

  % Finally optimize w.r.t h and sigma again
  nlmlF = @(t) normRotMargLikelihood( t(1)*oneVec, t(2)*oneVec, ...
    decomposition, A, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
  [~, optParams] = diRectWrap(nlmlF, diRectBounds, diRectOptions);
  sigmaSmOpts = optParams(1) * oneVec;
  sigmaPrOpts = optParams(2) * oneVec;

  % DEBUG
  fprintf('\t DEBUG: negMLF(A): %0.4f, negMLF(eye): %0.4f\n', negMLF(A), ...
    negMLF(eye(D,p)));

  % Finally Train the GP
  gpHyperParams.commonMeanFunc = commonMeanFunc;
  gpHyperParams.meanFuncs = meanFuncs;
  gpHyperParams.commonNoise = commonNoise;
  gpHyperParams.noises = noises;
  gpHyperParams.sigmaSms = sigmaSmOpts;
  gpHyperParams.sigmaPrs = sigmaPrOpts;
  % The function handles returned are for the transformed space but that is what
  % we want.
  Z = X * A; 
  [mu, KPost, Mus, KPosts, combinedZFuncH, funcHs] = ...
    addGPRegression(Z, y, Xtest, decomposition, gpHyperParams);
  combinedXFuncH = @(X) combinedZFuncH(X*A);

  % Finally return the learned decomposition
  learnedDecomp = cell(numGroups, 1);
  for i = 1:numGroups
    learnedDecomp{i} = permutedOrder(decomposition{i});
  end

end

