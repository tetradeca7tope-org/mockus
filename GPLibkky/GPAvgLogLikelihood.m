function [avgloglikl] = GPAvgLogLikelihood(mu, K, observation)
% Computes the average log-likelihood of an observation of a GP.

  n = size(mu,1);

  if isempty(K)
  % then return the negative sum of squared errors
    loglikl = - norm(mu - observation)^2;
    avgloglikl = loglikl / n;

  else
    maxdiag = max(diag(K));
    TOL_RATIO = 1e-3;
    % remove elements whose diagonals are too close to 0.
    idxs_to_keep = ones(n,1);
    for i = 1:n
      if K(i,i)/maxdiag < TOL_RATIO
        idxs_to_keep(i) = 0;
      end
    end
    num_discarded_pts = n - sum(idxs_to_keep);
    idxs_to_keep = idxs_to_keep > 0;
    mu = mu(idxs_to_keep);
    K = K(idxs_to_keep, idxs_to_keep);
    observation = observation(idxs_to_keep);
    num_pts = n - num_discarded_pts;

               
    loglikl = - (n/2)*log(2*pi) - 0.5*sum(log(diag(K))) + ...
              - (observation - mu)' * ( (1./ diag(K)) .* (observation - mu)) / 2;
              
    avgloglikl = loglikl/num_pts;
  end

end
