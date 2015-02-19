function A = numerical_2D_integration_wrap(th, z)

  N = sqrt(size(th,1));
  th1 = th(:,1); th2 = th(:,2);

  Th1 = reshape(th1, N, N);
  Th2 = reshape(th2, N, N);

  Z = reshape(z, N, N);
  A = numerical_2D_integration(Z, Th1, Th2);
  
end
