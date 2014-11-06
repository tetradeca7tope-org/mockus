% Unit test for truncated Gaussian Mean

muVals    = [0 2 4 10 200 -600]';
sigmaVals = [1 5 10 3 40  400]';
truncVals = [0 0 10 -8 230 -100]';

for i = 1:6
  mu = muVals(i);
  sigma = sigmaVals(i);
  trunc = truncVals(i);
  truncMean = truncatedGaussianMean(mu, sigma, trunc);

  % Now compute the qty via integration
  intFunc = @(t) (t  - trunc).* normpdf(t, mu, sigma);
  approxTruncMean = integral(intFunc, trunc, mu + 6*sigma);

  fprintf('Computed = %f\tApproximated = %f\n', truncMean, approxTruncMean);

end

% Check vectorized version
% truncVals = 0;
truncatedGaussianMean(muVals, sigmaVals, truncVals),
