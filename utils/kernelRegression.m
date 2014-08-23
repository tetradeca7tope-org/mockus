function fhat = kernelRegression(x, X, Y, h, k)
% Performs Kernel regression on 1D data.
% x: the point at which the kernel should be evaluated. X: data. Y: labels
% h: bandwidth, k: order of polynomial

  n = size(X, 1);

  diffs = X - x;
  omega = 1/sqrt(2*pi*h^2) * exp( - diffs.^2 / (2*h^2) );  
  Omega = diag(omega);

  U = zeros(n, k+1);
  U(:,1) = ones(n,1);
  for i = 2:(k+1)
    U(:,i) = diffs.^k;
  end
  
  theta_hat = pinv(U' * Omega * U) * (U' * Omega * Y);
  fhat = theta_hat(1);

end
