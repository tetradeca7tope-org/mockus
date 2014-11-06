function util = getUCBUtility(x, funcH, numEvals)

  % Prelims
  numDims = size(x, 2); % Expecting each x to be a row vector here.

  % Set beta_t. Using Recommendation from Section 6 in Srinivas et al, ICML 2010
  delta = 0.01; % something we need to set.
  t = numEvals + 1;
  beta_t = 2 * log( numDims * (t*pi)^2 / (6 * delta) ) / 5;

  % Obtain mean and standard deviation
  [mu, sigma] = funcH(x);
  util = mu + sqrt(beta_t) * sigma;

end

