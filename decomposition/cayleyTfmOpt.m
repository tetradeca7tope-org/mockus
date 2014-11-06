function [XOpt, optVal] = cayleyTfmOpt(func, Xinit, initVal, params)
% This function optimizes the function func subject to the constraint that
% X'X = Xinit'Xinit. (Usually, we init with Xinit'Xinit = I).
% func : A function handle which returns the to optimize compute the value of the objective.
%   If you are not using lineSearch this is not necessary.
% gradientF : A function handle that returns the gradient of func at X
% Xinit: Init point for the descent procedure.
% params: ancillary params for the problem

  % prelims
  n = size(Xinit, 1);
  p = size(Xinit, 2);
  maxNumTauViols = 5;

  if ~exist('params', 'var')
    params = struct;
  end

  if isfield(params, 'numIters'), numIters = params.numIters;
  else, numIters = 20;
  end

  if isfield(params, 'useLineSearch'), useLineSearch = params.useLineSearch;
  else, useLineSearch = true;
  end

  if ~useLineSearch
    if isfield (params, 'stepSize'), stepSize = params.stepSize;
    else stepSize = @(arg) 0.5/(arg);
    end
  end

  if isempty(initVal) & useLineSearch, initVal = func(Xinit);
  else initVal = [];
  end

  if ~isfield(params, 'initTau') | isempty(params.initTau)
    params.initTau = 1;
  end

  if ~isfield(params, 'XTol') | isempty(params.XTol)
    params.XTol = 0;
  end

  if isempty(Xinit)
  % Initialize to a random orthogonal matrix
    Xinit = rand(n, p); Xinit = orth(Xinit);
  end

  % Now start the iterative process
  X = Xinit; 
  val = initVal;
  tau = params.initTau;
  numTauViols = 0;

  % Choose best
  optVal = inf;

  for iter = 1:numIters
    
    % First compute the gradient
    [~, G] = func(X);
    A = G * X' - X * G';
%     PX = ( eye(n) - 0.5 * X * X'); PXG = PX * G; A = PXG * X' - X *PXG';
    if useLineSearch
      [X, val, tau] = lineSearch(func, A, X, val, max(tau, 2*params.XTol));
%       val, tau,
    else
      tau = stepSize(iter);
      X = getPointInDescentDirection(X, A, tau);
      val = func(X);
    end

    % Save the minimum
    if val < optVal
      optVal = val;
      XOpt = X;
    end

    % Quit if the change is too little
    if tau < params.XTol
      numTauViols = numTauViols + 1;
    end
    if numTauViols == maxNumTauViols, break;
    end

  end

end


% Uses back/forward tracking to do a line search
function [Y, val, tau] = lineSearch(func, A, Xcurr, currVal, initTau)

  n = size(Xcurr, 1);
  p = size(Xcurr, 2);

  % First step with tau = 1;
  tau = initTau;
  Y = Xcurr;
  Y = getPointInDescentDirection(Y, A, tau);
  val = func(Y);

  % Now back/ forward track
  if currVal < val
    while currVal < val
      tau = tau/2;
      Y = getPointInDescentDirection(Y, A, tau);
      val = func(Y);
      fprintf('tau: %0.4f val: %0.4f, currVAl: %0.4f\n', tau, val, currVal);
    end
  else
    while currVal > val
      tau = 2*tau;
      Y = getPointInDescentDirection(Y, A, tau);
      val = func(Y);
      fprintf('tau: %0.4f val: %0.4f, currVAl: %0.4f\n', tau, val, currVal);
    end
    tau = tau/2;
    Y = getPointInDescentDirection(Y, A, tau);
  end

end


% Obtains a point in the descent direction
function [Y] = getPointInDescentDirection(X, A, tau)
  n = size(A, 1);
  S = (eye(n) + tau * A/2);
  T = (eye(n) - tau * A/2);
%   S, T, A, tau, pause,
  Q = S \ T; % TODO: implement SMW here.
  Y = Q * X;
end

