% Unit test for bauesOptGP.m

close all;
clear all;
addpath ~/libs/kky-matlab/utils/
addpath ~/libs/kky-matlab/ancillary/
addpath ~/libs/kky-matlab/GPLibkky/

% choose function
funcChoice = 2;

% Function 1
if funcChoice == 1
  numDims = 1;
  func = @(t) t.^2 /4 + sin(7*t);
  bounds = [0 1];
  maxNumFunEvals = 15;

% Function 2
elseif funcChoice == 2
  numDims = 4;
  [func, funcProperties] = get3ModalFunction(numDims);
  bounds = funcProperties.bounds;
  maxNumFunEvals = min(8 * 2^numDims, 1000);
  maxPt = funcProperties.maxPt;
  fprintf('True maxVal, maxPt: %f, %s\n', func(maxPt), mat2str(maxPt) );

% Function 3
elseif funcChoice == 3
  numDims = 2;
  maxNumFunEvals = 30;
  bounds = [-2 2;-2 2];
  func = @(x) ...
    - (1+(x(:,1)+x(:,2)+1).^2.*(19-14.*x(:,1)+3.*x(:,1).^2 ...
        -14.*x(:,2)+6.*x(:,1).*x(:,2)+3.*x(:,2).^2))...
    .*(30+(2.*x(:,1)-3.*x(:,2)).^2.*(18-32.*x(:,1)+12.*x(:,1).^2 ...
        +48.*x(:,2)-36.*x(:,1).*x(:,2)+27.*x(:,2).^2));
end

  % 3.1 Call BO
  boParams.utilityFunc = 'EI';
  boParams.optPtStdThreshold = 0.002;
  boParams.alBWLB = 1e-4;
  boParams.alBWUB = 50;
  boParams.numInitPts = 4*numDims;
%   boParams.diRectParams.maxits = 8;
  boParams.gpNoiseLevel = 0.2; % * std(f(th));
  [maxVal, maxPt, boQueries, boVals, boHist] = ...
    bayesOptGP(func, bounds, maxNumFunEvals, boParams);

  % 3.2 Call DIRECT
  diRectOpts.maxevals = maxNumFunEvals;
  [~, ~, diRectHist] = diRectWrap(func, bounds, diRectOpts);
  numFuncEvaluations = diRectHist(end, 2);

  % 4. Plot iteration statistics
  figure;
  plot(boHist, 'b'); hold on,
  plot(diRectHist(:,2),diRectHist(:,3), 'r');
  legend('BO', 'DiRect');

  % 5. Plot the function
  if numDims == 1
    figure;
    res = 200;
    th = linspace(0,1, res);
    % Plot
    plot(th, func(th), 'b'); hold on,
    plot(boQueries, boVals, 'kx');
    plot(maxPt, maxVal, 'r*', 'MarkerSize', 10);
  elseif numDims == 2
    figure;
    plot2DFunction(@(t) (func(t)), bounds, 'contour'); hold on,
    plot(maxPt(1), maxPt(2), 'rx', 'MarkerSize', 10, 'LineWidth', 3);
    title('func');
  end

