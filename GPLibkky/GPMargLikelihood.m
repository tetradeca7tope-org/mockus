function [mu, K, funcH, sigmaSmOpt, sigmaPrOpt] = GPMargLikelihood(X, y, ...
  Xtest, hyperParams)
% Picks the hyper parameters for the GP using the marginal likelihood criterion.
% Structure for hyper parameters:
% - if hyperParams.sigmaSm is not 0, that will be the value used for the
%   smoothness. If it is zero then we do one of the following:
%    - If hyperParams.sigmaSmRange is a vector of length two we will optimize
%      between (sigmSmRange(1), sigmaSmRange(2)).
%    - If it is a scalar we will
%      optimize in the range logspace(-2, 2)*sigmaSmRange.
%    - If sigmaSmRange is empty, then we will set sigmaSmRange to the standard
%      deviation of x and adopt step 2.
% - Same applies to hyperParams.{sigmaPr, sigmaPrRange}. The default value for
%   sigmaPrRange would be the standard deviation of the y values
% - hyperParams.noise contains the observation noise. If not given we will add
%   some noise for numerical stability. Zero observation noise should be
%   specified.

  % prelims
  numDims = size(X, 2);
  numPts = size(X, 1);
  diRectOpt.maxits = 8; % max iterations for DiRect

  % Prepare bounds for the parameters
  % Bandwidth
  % ---------------------------------------------------
  if ~isempty(hyperParams.sigmaSmRange) 
    sigmaSmRange = hyperParams.sigmaSmRange;
  else
    stdX = norm(std(X, 0, 1))/sqrt(numDims);
    if stdX == 0, stdX = 1, end;
    sigmaSmRange = stdX;
  end
  % Scale 
  % ---------------------------------------------------
  if ~isempty(hyperParams.sigmaPrRange) 
    sigmaPrRange = hyperParams.sigmaPrRange;
  else
    sigmaPrRange = std(y);
  end
  % Now set the Bounds
  % ---------------------------------------------------
  if numel(sigmaSmRange) == 2, sigmaSmBounds = sigmaSmRange;
  else, sigmaSmBounds = [0.01 100] * sigmaSmRange;
  end
  if numel(sigmaPrRange) == 2, sigmaPrBounds = sigmaPrRange;
  else, sigmaPrBounds = [0.01 100] * sigmaPrRange;
  end

  % if no noise is given, add a little bit of noise
  if ~isfield(hyperParams, 'noise')
    hyperParams.noise = 0.01 * min(sigmaPrRange(1), std(y)) * ones(numPts, 1);
  end

  if hyperParams.sigmaSm ~= 0 && hyperParams.sigmaPr ~= 0
    sigmaSmOpt = hyperParams.sigmaSm;
    sigmaPrOpt = hyperParams.sigmaPr;

  elseif hyperParams.sigmaPr ~= 0 % only optimize over sigmaSm
    sigmaPrOpt = hyperParams.sigmaPr;
    % Prepare the function to be maximized by Direct
    nlmlF = @(t) normalizedMargLikelihood(t, sigmaPrOpt, X, y, ...
                                          hyperParams.noise);
    [~, sigmaSmOpt] = diRectWrap(nlmlF, sigmaSmBounds, diRectOpt);
    fprintf('Picked bw/scale = %.5f (%0.4f, %.4f), %.5f (fixed)\n', ...
                sigmaSmOpt, sigmaSmBounds(1), sigmaSmBounds(2), ...
                sigmaPrOpt );

  elseif hyperParams.sigmaSm ~= 0 % only optimize over sigmaPr
    sigmaSmOpt = hyperParams.sigmaSm;
    % Prepare the function to be maximized by Direct
    nlmlF = @(t) normalizedMargLikelihood(sigmaSmOpt, t, X, y, ...
                                          hyperParams.noise);
    [~, sigmaPrOpt] = diRectWrap(nlmlF, sigmaPrBounds, diRectOpt);
    fprintf('Picked bw/scale = %.5f (given), %.5f (%0.4f, %.4f)\n', ...
                sigmaSmOpt, ...
                sigmaPrOpt, sigmaPrBounds(1), sigmaPrBounds(2) );
    
  else  % optimize over both
    % Prepare the function to be maximized by Direct
    bounds = [sigmaSmBounds; sigmaPrBounds];
    nlmlF = @(t) normalizedMargLikelihood(t(1), t(2), X, y, hyperParams.noise);
    [~, optParams] = diRectWrap(nlmlF, bounds, diRectOpt);
    sigmaSmOpt = (optParams(1));
    sigmaPrOpt = (optParams(2));
    fprintf('Picked bw/scale = %.5f (%0.4f, %.4f), %.5f (%0.4f, %.4f)\n', ...
                sigmaSmOpt, sigmaSmBounds(1), sigmaSmBounds(2), ...
                sigmaPrOpt, sigmaPrBounds(1), sigmaPrBounds(2) );
  end

  % Finally retrain the GP 
  hyperParams.sigmaSm = sigmaSmOpt;
  hyperParams.sigmaPr = sigmaPrOpt;
  [mu, ~, K, funcH] = GPRegression(X, y, Xtest, hyperParams);

end

function nlml = normalizedMargLikelihood(sigmaSm, sigmaPr, X, y, noise)
  D = Dist2GP(X, X);
  K = sigmaPr * exp(-0.5*D/sigmaSm^2);
  Ky = K + diag(noise);
  nlml = -1/2 * ( y' * (Ky \ y) - log(det(Ky)) );
end

