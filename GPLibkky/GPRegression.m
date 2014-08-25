function [mu, stddev, Kpost, funcH] = GPRegression(X, y, Xtest, hyperParams, ...
                                              runtimeParams)
% Function for performing Gaussian Process Regression. Uses a Gaussian kernel.
% Updating this to use the Cholesky Decomposition on K11. The previous version
% is in GPRegressionOld.m

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
    % First prep the Kernel Matrix
    D11 = Dist2GP(X, X);
    K11 = sigmaPr * exp(-0.5*D11/sigmaSm^2);
    D12 = Dist2GP(X, Xtest);
    K12 = sigmaPr * exp(-0.5*D12/sigmaSm^2);
    D22 = Dist2GP(Xtest, Xtest);
    K22 = sigmaPr * exp(-0.5*D22/sigmaSm^2);

    % This routine does Cholesky Decomposition. If Cholesky fails due to
    % precision issues, then it adds an element to the diagonal.
    cholDecompSucc = false;
    K_ = K11 + diag(noise);
    diagPower = min( ceil(log10(abs(min(diag(K11)))))-1, -11);
    if ~(abs(diagPower) < inf)
      diagPower = -11;
    end
    % Now keep trying unitl Cholesky Decomp is successful
    while ~cholDecompSucc
      try
        L = chol(K_, 'lower');
        cholDecompSucc = true;
      catch err
        fprintf('CHOL failed in GP Regression with diagPower: %d\n', diagPower);
        diagPower = diagPower + 1;
        K_ = K_ + (10^diagPower)*eye(size(K11));
      end
    end
%     K_, (inv(L))'*inv(L), inv(L*L'),
    % Compute this value alpha too
    alpha = L' \ (L \ (y - meanX));
    
    % Obtain the outputs
    [mu, stddev, Kpost] = GPComputeOutputs(Xtest, X, L, alpha, ...
      sigmaPr, sigmaSm, meanFunc);

  else 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO: Better implementation for large amounts of training data
  end 

  if (runtimeParams.retFunc)
    funcH = @(Xte) ...
      GPComputeOutputs(Xte, X, L, alpha, sigmaPr, sigmaSm, meanFunc); 
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


function [yMu, yStd, yK] = GPComputeOutputs(Xte, Xtr, L, alpha, ...
  sigmaPr, sigmaSm, meanFunc)
% This function will be passed as a function handle to the output
% L is the Cholesky Decomp of the Kernel matrix.

  % Prelims
  meanXte = meanFunc(Xte);
  D12 = Dist2GP(Xtr, Xte);
  K12 = sigmaPr * exp(-0.5*D12/sigmaSm^2);
  D22 = Dist2GP(Xte, Xte);
  K22 = sigmaPr * exp(-0.5*D22/sigmaSm^2);

  % Now compute the outputs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 1. The predictive Mean
  yMu = meanXte + K12' * alpha;
  % 2. Predictive Variance
  V = L \ K12;
  yK = K22 - V'*V;
  yStd = sqrt(real(diag(yK)));
end

