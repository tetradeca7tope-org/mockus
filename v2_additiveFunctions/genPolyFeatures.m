function [Phi, derivs] = genPolyFeatures(X, genDerivs)

  if ~exist('genDerivs', 'var')
    genDerivs = false;
  end

  k = 2;

  n = size(X, 1);
  d = size(X, 2);

  Phi = [X X.^2 X.^3];

  combs = combnk(1:d, 2);
  N = size(combs, 1);
  newFeats = zeros(n, N);
  for j = 1:N
    newFeats(:,j) = prod( X(:, combs(j, :)), 2 );
  end
  Phi = [Phi newFeats];

%   for i = 2:k
%     combs = combnk( 1:d, i);
%     N = size(combs, 1);
%     newFeats = zeros(n, N);
%     for j = 1:N
%       newFeats(:,j) = prod( X(:, combs(j, :)), 2 );
%     end
%     Phi = [Phi newFeats];
%   end

  if genDerivs
    numTotalFeats = size(Phi, 2);
    derivs = zeros(numTotalFeats, d, n);
    for nIter = 1:n
      x = X(nIter, :)';
      currDerivs = zeros(numTotalFeats, d);
      currDerivs(1:d, :) = eye(d);
      currDerivs(d+1:2*d, :) = (2 * diag(x));
      currDerivs(2*d+1: 3*d, :) = (3 * diag(x.^2));

      combs = combnk(1:d, 2);
      N = size(combs, 1);
      for j = 1:N
        idx = 3*d + j;
        currDerivs(idx, combs(j, 1)) = x(combs(j, 2));
        currDerivs(idx, combs(j, 2)) = x(combs(j, 1));
      end

%     % Incomplete =============================
%     for i = 2:k
%       combs = combnk(1:d, i);
%       N = size(combs, 1);
%       for j = 1:N
%         feats(:,j) = prod( X(:, combs(j, :)), 2 );
%       end
%     end
%     % Incomplete =============================

      derivs(:,:,nIter) = currDerivs;
    end
  else
    derivs = [];
  end

end

