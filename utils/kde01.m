function [est_probs, f, h] = kde01(X, candidate_hs, doLogit)
% Performs Kernel density estimation on a dataset confined to [0,1]^d by first
% applying the logit transform. Returns the estimated probabilities and
% a function handle
% The candidate_hs need to be in the logit transformed R^d space

  if ~exist('candidate_hs', 'var')
    candidate_hs = [];
  end
  if ~exist('doLogit', 'var')
    doLogit = false;
  end

  if doLogit
    logit_X = logit(X);
    [~, flogit, h] = kde(logit_X, candidate_hs);
    f = @(data) flogit(logit(data)) ./ prod(data .* (1 - data), 2) ;
  else
    [~, f, h] = kde(X, candidate_hs);
  end
  % return estimated probs for the data
  est_probs = []; % f(X); doing this for now to speed things up
end
