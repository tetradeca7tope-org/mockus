function A = numerical_2D_integration(Z, X, Y)

  n1 = size(Z, 1);
  n2 = size(Z, 2);

  A = 0;
  for i = 1:(n1-1)
    for j = 1:(n2-1)
      accum = (1/4) * (X(1,j+1) - X(1,j)) * (Y(i+1,1) - Y(i,1)) * ...
              (Z(i,j) + Z(i+1,j) + Z(i,j+1) + Z(i+1,j+1));
      if abs(accum) < inf
        A = A + accum;
      end
    end
  end

end
