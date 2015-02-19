function [est_probs, f, h] = kde(X, candidate_hs)
% Performs Kernel density estimation on the data. Picks the bandiwidth using
% 10-fold cross validation. Assumes domain [-inf, inf]^d
% If domain is [0, 1]^d use kde01
% input - X: data in an nxd matrix. 
% output - est_probs: estimated probabilities, f: A function handle to estimate
% the density on new data. h: optimal bandwidth

  n = size(X, 1);
  d = size(X, 2);
  num_partitions_kfoldcv = min(10, n);

  % First determine the candidates for cross validation if they were not given
  if ~exist('candidate_hs', 'var') | isempty(candidate_hs)
    stdX = std(X);
    if stdX == 0, stdX = 1; end
    silverman_h = 1.06 * stdX / n^(-1/5);
    candidate_hs = logspace(-2,2,10)' * silverman_h;
  end
  num_cands = size(candidate_hs, 1);

  % Use cross validation to pick the right h
  if num_cands == 1 % if only one point given, use that
    best_h = candidate_hs;

  else % otherwise cross-validate

    % If only 1 point, choose the largest bandwidth 
    if n == 1 % if you have only one point
      best_h = candidate_hs(end);
    else
      % shuffle the data
      X = X(randperm(n), :);
      best_likl = -inf;
      best_h = candidate_hs(1);

      for cand_iter = 1:num_cands
        curr_h = candidate_hs(cand_iter);
        cv_loglikl = KFoldExperiment(X, curr_h, num_partitions_kfoldcv);
        if cv_loglikl >= best_likl
          best_likl = cv_loglikl;
          best_h = curr_h;
        end
      end % for cand_iter
    end % if n ==1

  end % num_cands==1

  % finally prep results for returning
  h = best_h;
  f = @(data) kdeEstAtPts(data, X, h); % ( sum(GaussKernel(h, X, data))'/ n); 
%   est_probs = f(X); % doing this to speed things up
  est_probs = []; 

end


function [cv_loglikl] = KFoldExperiment(X, h, num_partitions_kfoldcv)

  loglikl_accum = 0;
  n = size(X, 1);

  for kfold_iter = 1:num_partitions_kfoldcv
    % Set the partition up.
    test_start_idx = round( (kfold_iter-1)*n/num_partitions_kfoldcv + 1 );
    test_end_idx   = round( kfold_iter*n/num_partitions_kfoldcv );
    train_indices = [1:test_start_idx-1, test_end_idx+1:n]';
    test_indices = [test_start_idx : test_end_idx]';
    num_test_data = test_end_idx - test_start_idx + 1;
    num_train_data = n - num_test_data;
    Xtr = X(train_indices, :);
    Xte = X(test_indices, :);
    % Compute the log likelihood
    Pte = kdeEstAtPts(Xte, Xtr, h); % sum(GaussKernel(h, Xtr, Xte))' / n;
    avg_log_likl = sum(log(Pte)) / num_test_data;
    loglikl_accum = loglikl_accum + avg_log_likl;
  end
  cv_loglikl = loglikl_accum / num_partitions_kfoldcv;
end


function ests = kdeEstAtPts(pts, X, h)
  numPts = size(pts, 1);
  numData = size(X, 1);
  maxNumPts = max(1e7, numData);
  ptsPerPartition = min( numPts, ceil(maxNumPts/numData));

  % Now do the estimation
  cumNumPts = 0;
  ests = zeros(numPts, 1);
  while cumNumPts < numPts
    currNumPts = min(ptsPerPartition, numPts - cumNumPts);
    ests(cumNumPts+1: cumNumPts + currNumPts) = ...
      mean(GaussKernel(h, X, pts(cumNumPts+1: cumNumPts + currNumPts, :)))';
    cumNumPts = cumNumPts + currNumPts;
  end
end

% function est = kdeEstAtPt(pt, X, h)
%   est = mean(GaussKernel(h, X, pt));
% end
