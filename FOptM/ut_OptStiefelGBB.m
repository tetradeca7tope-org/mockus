function ut_OptStiefelGBB
% My unit test for OptStiefelGBB.m

%   n = 500; p = 300;
  n = 5; p = 3;

  B = eye(p,n);
  Ainit = rand(n, p); Xinit = orth(Ainit);

  obj = @(arg) func(arg, B);

  % Print out initial values
%   Xinit,
  initVal = obj(Xinit),

  opts.record = 1;
  opts.mxitr  = 1000;
  opts.xtol = 1e-5;
  opts.gtol = 1e-5;
  opts.ftol = 1e-8;
  out.tau = 1e-3;
  [X, out] = OptStiefelGBB(Xinit, obj, opts);

%   X,
  out,
  


end


function [F, G] = func(X, B)

  F = trace(X * B);
  G = B';

end
