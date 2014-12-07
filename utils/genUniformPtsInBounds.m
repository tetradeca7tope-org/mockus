function X = genUniformPtsInBounds(numPts, bounds)
% Generates points from the uniform distribution in the given bounds
% The first column is the lower bound on each dimension and the second is the
% upper bound.
  numDims = size(bounds, 1);
  U = rand(numPts, bounds);
  X = bsxfun(@plus, bsxfun(@times, U, bounds(:,2) - bounds(:,1)), bounds(:,1));
end

