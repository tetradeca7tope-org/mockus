function K = ...
    combinedKernelNoise(X, decomposition, bws, scales, noises, commonNoise)
  % Same as complete kernel but with noise
  K = combinedKernel(X, X, decomposition, bws, scales) + ...
    (sum(noises) + commonNoise) * eye(size(X, 1));
end
