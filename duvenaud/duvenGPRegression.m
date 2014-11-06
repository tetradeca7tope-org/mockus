function [mu, stddev, Kpost, funcH] = duvenGPRegression(X, y, Xtest, ...
  hyperParams, runtimeParams)
% This implements Additive GP Regression as described in Duvenaud et al,
% Additive Gaussian Processes, NIPS 2012

  % Prelims
  numTrData = size(X, 1);
  numDims = size(X, 2);

  
  % Set the hyper parameters for each GP
  % Mean Function
  % ---------------------------------------------------------------------------
  if isempty(hyperParams.meanFunc)
    meanFunc = @(arg) mean(y) * ones( size(arg, 1), 1);
  else
    meanFunc = hyperParams.meanFunc;
  end

  % noise parameters
  if isempty(hyperParams.noise)
    noise = zeros(numTrData, 1);
  elseif numel(hyperParams.noise) == 1
    noise = hyperParams.noise * ones(numTrData, 1);
  else
    noise = hyperParams.noise;
  end


end

