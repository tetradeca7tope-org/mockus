function [Phi, Derivs, groups, derivGroups] = ...
  genDecompPolyFeatures(X, A, d, genDerivs)
% This first applies the transformation AX and then generates polynomial
% features in groups of D

  Z = X * A;
  D = size(Z, 2);
  n = size(Z, 1);

  if ~exist('genDerivs', 'var')
    genDerivs = false;
  end

  % Initialize features
  Phi = ones(n, 1);
  Derivs = zeros(1, D, n);
  numGroups = ceil(D/d);
  phiGroups = {};
  derivGroups = {};

  for i = 1:numGroups
    start_idx = (i-1)*d + 1; pre = start_idx-1; %1:(start_idx-1),
    end_idx = min((d*i), D); post = D - end_idx; % (end_idx+1):D,
    Zcurr = Z(:, start_idx:end_idx);
    [currPolyFeats, currPolyDerivs] = genPolyFeatures(Zcurr, genDerivs);

    Phi = [Phi currPolyFeats];
    groups{i} = currPolyFeats;

    if genDerivs
      numPolyFeats = size(currPolyFeats, 2);
%       size(currPolyDerivs), pre, post, size(currPolyFeats),
      currPolyDerivs = [zeros(numPolyFeats, pre, n) currPolyDerivs, ...
        zeros(numPolyFeats, post, n)];
%       size(currPolyDerivs),
      Derivs = [Derivs; currPolyDerivs];
      derivGroups{i} = currPolyDerivs;
    else
      Derivs = [];
      derivGroups = [];
    end
  end

end

