function [betaOpt, AOpt] = optimizeABeta(Y, X, params)
% Optimizes for beta and A via a coordinate descent scheme.

  if ~isfield(params, 'numIters') | isempty(params.numIters), numIters = 10;
  else, numIters = params.numIters;
  end

  if ~isfield(params, 'numInits') | isempty(params.numInits), numInits = 10;
  else, numInits = params.numInits;
  end
  
  if ~isfield(params, 'verbose') | isempty(params.verbose), verbose = false;
  else, verbose = params.verbose;
  end

  if ~isfield(params, 'maxOptAIters') | isempty(params.maxOptAIters)
    maxOptAIters = 10;
  else
    maxOptAIters = params.maxOptAIters;
  end

  % prelims
  numData = size(X, 1);
  numDims = size(X, 2);
  d = params.d;
  lambda = params.lambda;
  
  % Params for optimizing A
  stiefelOpts.record = 0;
  stiefelOpts.xtol = 1e-5;
  stiefelOpts.gtol = 1e-5;
  stiefelOpts.ftol = 1e-8; 
  stiefelOpts.tau = 1e-3;

  currBest = inf;

  for initIter = 1:numInits
    % Initialize A
    A = randn(numDims); A = orth(A);
    progress = zeros(numIters, 1);

    for iter = 1:numIters

      % First generate the features and optimize for beta
      Phi = genDecompPolyFeatures(X, A, d, false);
      numAugDims = size(Phi, 2);
      
      % Now optimize for beta
      beta = (Phi' * Phi + lambda * eye(numAugDims)) \ Phi' * Y;
      
      % Now optimize for A
      if iter == numIters, stiefelOpts.mxitr = min(15*maxOptAIters, 100);
      else, stiefelOpts.mxitr = maxOptAIters;
      end
      obj = @(arg) objA(arg, X, beta, Y, d);
      tic,
      [A, out] = OptStiefelGBB(A, obj, stiefelOpts);
      toc,

      progress(iter) = norm(Phi*beta - Y);
    end

      plot(progress); pause;

    % Finally reoptimize for beta
      Phi = genDecompPolyFeatures(X, A, d, false);
      numAugDims = size(Phi, 2);
      beta = (Phi' * Phi + lambda * eye(numAugDims)) \ Phi' * Y;

    % Now check against best value
    currVal = norm(Phi*beta - Y);
      % Print some results
      if verbose, initIter, beta, A, currVal,
      end
    if currVal < currBest
      AOpt = A;
      betaOpt = beta;
      currBest = currVal;
    end

  end

  if verbose
    fprintf('Best Value: %f\n', currBest);
    AOpt, betaOpt,
  end
end



function [F, G] = objA(A, X, beta, Y, d)
% Returns the value of the function and the gradient at A

  % generate features
  [Phi, dPhiDZ]  = genDecompPolyFeatures(X, A, d, true);
%   Phi,
%   size(Phi), size(dPhiDZ),
  diff = Phi * beta - Y;
  F = norm(diff)^2;

  % Some notation
  N = size(Phi, 2);
  D = size(A, 1);
  p = size(A, 2);
  n = size(X, 1);
  beta3 = reshape(beta, 1, 1, numel(beta));

  % Now we need to compute the derivatives
  G = zeros(D, p);
  for nIter = 1:n
    dPsidZ = dPhiDZ(:,:,nIter)' * beta;
    x = X(nIter, :)';
    currDPsiDA = bsxfun(@times, repmat(x, 1, p), dPsidZ');
    G = G + 2 * diff(nIter) * currDPsiDA;
  end

end

%   Slower method for computing the derivative
%   % Now we need to compute the derivatives
%   G = zeros(D, p);
%   for nIter = 1:n
%     currDPhiDA = zeros(D, p, N); 
%     x = X(nIter, :)';
%     for i = 1:N
% %     N, nIter, i, size(dPhiDZ), size(x), 
%       currDPhiDA(:,:,i) = bsxfun(@times, repmat(x, 1, p), dPhiDZ(i, :, nIter));
%     end
%     currDPhiDABeta = sum( bsxfun(@times, currDPhiDA, beta3), 3); 
%     % Finally multiply by diff(nIter)
% %     DPhiDABeta = diff(nIter) * currDPhiDABeta;
%     G = G + 2 * diff(nIter) * currDPhiDABeta;
%   end

