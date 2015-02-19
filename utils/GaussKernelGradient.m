function G = GaussKernelGradient(h, X, Y)
% Returns the Gradient of the Gaussian Kernel with respect to the points in X.
% returns a Dx nY x nX tensor 

  if ~exist('Y', 'var') | isempty(Y)
    Y = X;
  end

  % prelims
  nX = size(X, 1);
  nY = size(Y, 1);
  D  = size(X, 2);

  K = GaussKernel(h, X, Y);

  G = zeros(D, nY, nX);
  for i = 1:nX
    diffs = bsxfun(@minus, Y', X(i,:)');
    G(:, :, i) = bsxfun( @times, diffs ,K(i, :)) / h^2;
  end

end

