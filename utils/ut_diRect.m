% A unit test for diRect.m
% DiRect (Dividing Rectangles) is a zeroth order Branch & Bound like Global
% Optimzation algorithm.

clear all;
close all;

% Test 1
%%%%%%%%
fprintf('Test 1\n===================\n');
  % 1. Establish bounds for variables
%   options.maxits = 8;
  options.maxevals = 100;

% choose function
funcChoice = 1;

% Function 1
if funcChoice == 1
  numDims = 5;
  bounds = repmat([0, 1], numDims, 1);
  funcProperties.gaussVar = 0.005 * sqrt(numDims);
  funcProperties.covarDiag = funcProperties.gaussVar * ones(1, numDims);
  funcProperties.centres12 = [0.81 0.31; 0.59 0.79; 0.21 0.22];
  funcProperties.centresRest = [0.17 0.41 0.81]';
  funcProperties.mixProportions = [0.2 0.5 0.8]';
  funcProperties.centres = ...
    [funcProperties.centres12, repmat(funcProperties.centresRest, 1, numDims -2)];
  func = @(t) log( ...
    0.2 * mvnpdf(t, funcProperties.centres(1, :), funcProperties.covarDiag ) + ...
    0.5 * mvnpdf(t, funcProperties.centres(2, :), funcProperties.covarDiag ) + ...
    0.3 * mvnpdf(t, funcProperties.centres(3, :), funcProperties.covarDiag ) );
  maxPt = funcProperties.centres(2, :);
  fprintf('True maxVal, maxPt: %f, %s\n', func(maxPt), mat2str(maxPt) );

elseif funcChoice == 2
  numDims = 2;
  bounds = [-2 2;-2 2];
  func = @(x) ...
    - (1+(x(:,1)+x(:,2)+1).^2.*(19-14.*x(:,1)+3.*x(:,1).^2 ...
        -14.*x(:,2)+6.*x(:,1).*x(:,2)+3.*x(:,2).^2))...
    .*(30+(2.*x(:,1)-3.*x(:,2)).^2.*(18-32.*x(:,1)+12.*x(:,1).^2 ...
        +48.*x(:,2)-36.*x(:,1).*x(:,2)+27.*x(:,2).^2));
end

  % 3. Call DIRECT
  tic,
  [fmax, xmax, hist] = diRectWrap(func, bounds, options);
  toc,
  numFuncEvaluations = hist(end, 2);
  fprintf('fmax: %f, xmax: %s\n# func evals: %d\n', fmax, mat2str(xmax), ...
    numFuncEvaluations);

  % 4. Plot iteration statistics
  plot(hist(:,2),hist(:,3))
  xlabel('Fcn Evals');
  ylabel('f_{min}');
  title('Iteration Statistics for GP test Function');

  % Also plot the function out
  if numDims == 2
  figure;
    plot2DFunction(@(t) (func(t)), bounds, 'contour'); hold on,
    plot(xmax(1), xmax(2), 'rx', 'MarkerSize', 10, 'LineWidth', 3);
    title('log f');
  end

  % Randomly sample function evaluations and report the maximum
  Xrand = rand(10*numFuncEvaluations, numDims);
  tic,
  Yrand = func(Xrand);
  [maxVal, maxIdx] = max(Yrand);
  toc,
  [minVal, minIdx] = min(Yrand);
  fprintf('Rand maxVal, maxPt: %f, %s\n', maxVal, mat2str(Xrand(maxIdx, :)) );
  fprintf('minVal, minPt: %f, %s\n', minVal, mat2str(Xrand(minIdx, :)) );


  % 3. Call DIRECT

