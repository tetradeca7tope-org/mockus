function ut_cayleyTfmOpt
% My unit test for OptStiefelGBB.m

%   n = 500; p = 300;
  n = 5; p = 5;

  B = eye(p,n);
  Ainit = rand(n, p); Xinit = orth(Ainit);

  obj = @(arg) func(arg, B);
%   obj = @func2;
%   obj = @func3;

  % Print out initial values
%   Xinit,
  initVal = obj(Xinit),

  opts.record = 1;
  opts.mxitr  = 1000;
  opts.xtol = 1e-5;
  opts.gtol = 1e-5;
  opts.ftol = 1e-8;
  out.tau = 1e-3;
  tic,
  [X, out] = OptStiefelGBB(Xinit, obj, opts);
  toc,

  fprintf('Stiefel Implementation\n');
  X, out.fval,

  params.numIters = 10;
  params.XTol = 1e-5;
  params.initTau = 10;
  tic,
%   [X, optval] = cayleyTfmOpt(obj, Xinit, [], params);
%   toc,
% 
%   fprintf('My Implementation\n');
%   X, optval,
  
  obj(eye(n,p)),


end


function [F, G] = func(X, B)
  F =  - trace(X * B);
  G = - B';
end


function [F, G] = func2(X)
  F = norm(X, 'fro'); 
  G = 2*X;
end

function [F, G] = func3(X)
  F = norm(X, 'fro')^2; 
  G = 2*X * sqrt(F);
end
