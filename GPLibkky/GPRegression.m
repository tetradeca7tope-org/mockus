function [mu, stddev, K, funcH] = GPRegression(X, y, Xtest, hyperParams, ...
                                              runtimeParams)
% Function for performing Gaussian Process Regression. Uses a Gaussian kernel.

  sigmaSm = hyperParams.sigmaSm;
  sigmaPr = hyperParams.sigmaPr;
  meanFunc = hyperParams.meanFunc;
  noise = hyperParams.noise;
  num_tr_data = size(X, 1);
  num_te_data = size(Xtest, 1);

  % Run time parameters
  if ~exist('runtimeParams', 'var'), runtimeParams = struct(); end
  if ~isfield(runtimeParams, 'plotOn'), runtimeParams.plotOn = false; end
  if ~isfield(runtimeParams, 'retFunc'), runtimeParams.retFunc = false; end

  if isempty(noise)
    noise = zeros(size(X,1), 1);
  end
  if isempty(meanFunc) % if meanFunc is empty use the mean of the y points
    meanFunc = @(arg) mean(y);
  end

  % Precompute some preliminaries
  meanX = meanFunc(X);
  mean_Xte = meanFunc(Xtest);
  if size(meanX, 1) == 1
    meanX = meanX * ones(num_tr_data, 1);
    mean_Xte = mean_Xte * ones(num_te_data, 1);
  end
  if num_tr_data < inf
%     fprintf('Inverting full matrix !\n');
    D11 = Dist2GP(X, X);
    K11 = sigmaPr * exp(-0.5*D11/sigmaSm^2);
    D12 = Dist2GP(X, Xtest);
    K12 = sigmaPr * exp(-0.5*D12/sigmaSm^2);
    D22 = Dist2GP(Xtest, Xtest);
    K22 = sigmaPr * exp(-0.5*D22/sigmaSm^2);

    % obtain outputs
    % concatenate these 2 matrices so that we invert the matrix only once.
    QQ = (K11 + diag(noise)) \ [K12 (y - meanX)];
    mu = mean_Xte + K12' * QQ(:, end);
    
    K = K22 - K12' * QQ(:,1:end-1);
    stddev = real(sqrt(diag(K)));

  else  % num_data < 100
%     fprintf('Inverting 100x100 matrix!\n');
    % If num_tr_data is too large, then just do the GP on the closest 100 pts.
    % compute the predictions for each point separately.
    mu = zeros(num_te_data, 1);
    stddev = zeros(num_te_data, 1);
    K = []; % don't return this. difficult ot compute ?

    for te_iter = 1:num_te_data
      xte = Xtest(te_iter, :)';
      dist2_to_tr = Dist2GP(X, xte'); 
      [~, sorted_idxs] = sort(dist2_to_tr);
      nbd = sorted_idxs(1:100);
      xtr = X(nbd, :);

      % now do what you did before on xtr
      d11 = Dist2GP(xtr, xtr);
      k11 = sigmaPr * exp(-0.5*d11/sigmaSm^2);
      d12 = dist2_to_tr(nbd);
      k12 = sigmaPr * exp(-0.5*d12/sigmaSm^2);

      % obtain output
      qq = (k11 + diag(noise(nbd))) \ [k12 (y(nbd) - meanX(nbd))];
      mu(te_iter) = mean_Xte(te_iter) + k12' * qq(:, end);

      % return the covariance if necessary
      d22 = 0;
      k22 = sigmaPr;
      stddev(te_iter) = k22 - k12' * qq(:, 1:end-1);
    end % for te_iter

  end % else num_data < 100

  if (runtimeParams.retFunc)
  % Create a function handle so that the matrix need not be inverted all the
  % time
    Kmat = K11 + diag(noise);
    if rcond(Kmat) < 1e-16
      invMat = inv(Kmat);
    else
      invMat = pinv(Kmat);
    end
    funcH = @(Xte) GPFuncHandle(Xte, X, y, invMat, sigmaPr, ...
                                sigmaSm, meanFunc); 
  else
    funcH = [];
  end

  if (runtimeParams.plotOn && (size(X,2) ==1) )
    plot(X, y, 'kx', 'MarkerSize', 10); hold on,
    plot(Xtest, mu);
    plot(Xtest, mu + 3*stddev, 'g--');
    plot(Xtest, mu - 3*stddev, 'g--');
    title('GP results: original points(kx), estimated vals(b),error(g--)');
  end

end

% This function will be passed as a function handle to the output
function [yte, ystd, yK] = GPFuncHandle(Xte, Xtr, Ytr, invK, sigmaPr, sigmaSm, meanFunc)
  meanXte = meanFunc(Xte);
  meanXtr = meanFunc(Xtr);

  D12 = Dist2GP(Xtr, Xte);
  K12 = sigmaPr * exp(-0.5*D12/sigmaSm^2);
  D22 = Dist2GP(Xte, Xte);
  K22 = sigmaPr * exp(-0.5*D22/sigmaSm^2);
  yte = meanXte + K12' * invK * (Ytr - meanXtr);
  yK = K22 - K12' * invK * K12;
  ystd = real(sqrt(diag(yK)));
end
