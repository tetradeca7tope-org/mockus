function [mu, K, sigmaSmOpt, sigmaPrOpt] = GPKFoldCV(X, y, Xtest, ...
  num_partitions, candidates, hyperParams)
% Performs K-fold CV for the GP. candidates is a struct which gives the
% possible values for the smoothness and scale parameters of the GP.
% The number of partitions (K for K-fold CV) is given in the arg num_partitions
% If num_partitions is empty use K = min(num_data, 10);

  if isstruct(candidates)
    n1 = size(candidates.sigmaSmVals, 1);
    n2 = size(candidates.sigmaPrVals, 1);
    a1 = repmat(candidates.sigmaSmVals, n2, 1);
    a2 = repmat(candidates.sigmaPrVals, 1, n1)'; a2 = a2(:);
    cand_vals = [a1 a2];
%     cand_vals,
  else
    cand_vals = candidates;
  end

  if size(cand_vals, 1) == 1
    sigmaSmOpt = cand_vals(1);
    sigmaPrOpt = cand_vals(2);
  else
    % if the Cost function is not specified, the set it to be sum of squares
    if ~isfield(hyperParams, 'costFunc')
      hyperParams.costFunc = @(y1, y2) (y1 - y2).^2;
    end

    if isempty(num_partitions)
      num_partitions = min(size(X,1), 10);
    end
    % if fewer partitions than data fix that
    num_partitions = min(num_partitions, size(X,1));

    % Peform K-fold cross Validation
    num_cand_combinations = size(cand_vals);
    best_likl = -inf;
    for cand_iter = 1:num_cand_combinations
      hyperParams.sigmaSm = cand_vals(cand_iter, 1);
      hyperParams.sigmaPr = cand_vals(cand_iter, 2);
      curr_likl = KFoldExperiment(X, y, num_partitions, hyperParams);
      if curr_likl >= best_likl
        best_likl = curr_likl;
        sigmaSmOpt = cand_vals(cand_iter, 1);
        sigmaPrOpt = cand_vals(cand_iter, 2);
      end
    end
  end

  % Finally rerun GP for the entire test set.
  hyperParams.sigmaSm = sigmaSmOpt;
  hyperParams.sigmaPr = sigmaPrOpt;
  [mu, ~, K] = GPRegression(X, y, Xtest, hyperParams);

end


function [kfold_likl] = KFoldExperiment(X, y, num_partitions, hyperParams)
%% This function performs K-fold cross validation for the current set of
%% values and returns the average likelihood.
  
  m = size(X, 1);
  likl_accum = 0;
  % Shuffle the data
  data_ordering = randperm(m);
  X = X(data_ordering, :);
  y = y(data_ordering, :); 
  orig_noise = hyperParams.noise(data_ordering, :);

  for kfold_iter = 1:num_partitions
    test_start_idx = round( (kfold_iter-1)*m/num_partitions + 1 );
    test_end_idx   = round( kfold_iter*m/num_partitions );
    train_indices = [1:test_start_idx-1, test_end_idx+1:m];
    test_indices = [test_start_idx : test_end_idx];
    Xtr = X(train_indices, :);
    ytr = y(train_indices);
    Xte = X(test_indices,:);
    yte = y(test_indices);
    hyperParams.noise = orig_noise(train_indices);
    [mu, ~, K] = GPRegression(Xtr, ytr, Xte, hyperParams);
%     curr_avg_log_likl = GPAvgLogLikelihood(mu, K, yte);
    curr_avg_log_likl = - sum(hyperParams.costFunc(mu, yte));
    if isnan(curr_avg_log_likl)
      curr_avg_log_likl = -inf;
    end
    likl_accum = likl_accum + curr_avg_log_likl;
  end
  kfold_likl = likl_accum / m;
end

