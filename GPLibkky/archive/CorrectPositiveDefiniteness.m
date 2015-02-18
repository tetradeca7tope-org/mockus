function A = CorrectPositiveDefiniteness(A)

  [V, D] = eig(A);
  d = real(diag(D));
  V = real(V);
  d = d .* double(d > 0); 
  min(d),
  A = V * diag(d) * V';

end
