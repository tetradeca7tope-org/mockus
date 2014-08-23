function [l2] = l2divergence(X, Y)
% Estimates the L2 divergence between X & Y.
% TODO: kde creates a nxn matrix so don't use this for a large number of pts.

  n = size(X, 1);
  m = size(Y, 1);

  % shuffle the data
  X = X(randperm(n), :);
  Y = Y(randperm(m), :);

%   SPLIT_DATA = false;
  SPLIT_DATA = true;
  % TODO: need to data split. Otherwise, it doesn't work. Why ?
  if SPLIT_DATA
    n1 = round(n/2);
    n2 = n - n1;
    m1 = round(m/2);
    m2 = m - m1;
    X1 = X(1:n1, :);
    X2 = X(n1+1:end, :);
    Y1 = Y(1:m1, :);
    Y2 = Y(m1+1:end, :);

  else
    n1 = n;
    n2 = n;
    m1 = m;
    m2 = m;
    X1 = X;
    X2 = X;
    Y1 = Y;
    Y2 = Y;
  end

  % Obtain KDE estimates
  h = 0.01;
  [~, phat, hX] = kde(X1, h); 
  [~, qhat, hY] = kde(Y1, h); 

  est_int_pp = mean(phat(X2));
  est_int_qq = mean(qhat(Y2));
  est_int_pq = mean(phat(Y2));
  est_int_qp = mean(qhat(X2));

  l2 = sqrt( est_int_pp + est_int_qq - est_int_pq - est_int_qp);

end
