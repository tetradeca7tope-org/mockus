function [A, queries, vals, samples] = ccaDecomposition(func, numDims, ...
  numSamples, mcmcPropStd)
% Uses Barnabas' Cost Components Analysis method to identify the decomposition.

  % First use MCMC to obtain samples
  mcmcInitPt = 0.1 * randn(numDims, 1);

  % Now obtain the samples
  [samples, queries, vals] = ...
    customMCMC(numSamples, mcmcPropStd, mcmcInitPt, func);
%   numel(unique(samples)),

%   size(samples),
  % Now do ICA
  [AT, ~] = fastica(samples');
  A = AT';
  colNorms = sqrt( sum( A.^2) );
  A = bsxfun(@rdivide, A, colNorms);
  
end

