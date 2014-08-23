function V = GPDrawSamples(mu, K, numSamples)
% A function for drawing samples from a Gaussian Process
% V is a num_samples x m matrix, each row corresponding to a sample

  num_pts = size(mu,1);
  V = zeros(numSamples, num_pts);

  chol_decom_successful = false;
  diag_power = min(ceil( log10(abs(min(diag(K)))) ) - 1, -10);
  if ~(abs(diag_power) < inf)
    diag_power = -10;
  end

  while ~ chol_decom_successful
    try
      K_ = K + (10^diag_power)*eye(size(K));
      diag_power = diag_power + 1;
      chK = chol(K_);
      chol_decom_successful = true;
    catch err
      fprintf('CHOL failed with diag_power: %d\n', diag_power);
    end
  end

  U = randn(num_pts, numSamples);
  V = real(bsxfun(@plus, chK'*U, mu));
  V = V'; 

end
