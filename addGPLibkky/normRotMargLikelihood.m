function [nlml, AG] = normRotMargLikelihood(sigmaSms, sigmaPrs, decomposition, A, ...
  X, y, meanFuncs, commonMeanFunc, noises, commonNoise)
% Returns the normalized marginal likelihood and if requested its derivative
% w.r.t A.
% Decomposition is the decomposition after applying Z = X*A;

  % prelims
  numPts = size(X, 1);
  numGroups = numel(decomposition);
  D = size(A, 1);
  p = size(A, 2);

  % apply transformation and compute normalized marginal likelihood
  Z = X * A;
%   Z, decomposition, sigmaSms, sigmaPrs, noises, commonNoise,
  Ky = combinedKernelNoise(Z, decomposition, sigmaSms, sigmaPrs, noises, ...
        commonNoise);
  L = stableCholesky(Ky);
  y_ = y - combinedRotMeanFunc(X, Z, commonMeanFunc, meanFuncs, decomposition);
  alpha = L' \ (L \ y);
  nlml = -1/2 * y_' * alpha - sum(log(diag(L))) - numPts/2 * log(2*pi);

  if nargout == 2
  % Then also compute the gradient w.r.t. A
    AG = zeros(D, p);

      for k = 1:numGroups

        % First compute the kernel for this group
        coords = decomposition{k};
        sm = sigmaSms(k);
        pr = sigmaPrs(k);
        Kk = augKernel(Z, Z, coords, sm, pr);
        Zk = Z(:, coords);

        % Now go through each of the coordinates in this group
        for j = 1:numel(decomposition{k})
          currCoord = decomposition{k}(j);
          dKdz1 = GaussKernelGradientCoords(sm, Zk, Zk, Kk, j);

          parfor i = 1:D

            % First construct the nxn matrix dK/dAij
            diffs = bsxfun(@minus, repmat(X(:,i), 1, numPts), X(:,i)');
            dKdAij = dKdz1 .* diffs;

            % Now compute 1/2 * tr[ (alpha*alpha' - inv(Ky)) dK/dAij ]
            t1 = (alpha' * dKdAij) * alpha;
            t2 = trace( L'\ (L \ dKdAij) );
            AG(i, currCoord) = 0.5 * (t1 - t2);
          
          end
        end
      end
  end

end


function mu0 = combinedRotMeanFunc(X, Z, commonMeanFunc, meanFuncs, ...
                 decomposition)
% The common Mean Func takes in X as its arguments while meanFuncs take in Z =
% X*A as its arguments

  mu0 = commonMeanFunc(X);
  numGroups = numel(decomposition);
  for k = 1:numGroups
    coords = decomposition{k};
    mu0 = mu0 + meanFuncs{k}( Z(:, coords) );
  end

end
