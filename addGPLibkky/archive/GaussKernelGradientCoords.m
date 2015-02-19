function G = GaussKernelGradientCoords(bw, X, Y, K, coords)
  % computes the gradient with respect to the given coords
  % K is the Gaussian Kernel Matrix

  % prelims
  nX = size(X, 1);
  nY = size(Y, 1);
  D  = size(X, 2);
  if isempty(coords), coords = (1:D)';
  end
  numCoords = numel(coords);

  G = zeros(nX, nY, numCoords);
  for i = 1:numCoords
    coordDiffs = bsxfun(@minus, repmat(X(:,i), 1, nY), Y(:,i)');
    G(:,:,i) = - K .*coordDiffs / bw^2;
  end

end

