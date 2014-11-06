% Unit test for bauesOptGP.m

close all;
clear all;
addpath ~/libs/kky-matlab/utils/
addpath ~/libs/kky-matlab/ancillary/
addpath ~/libs/kky-matlab/GPLibkky/
addpath ../BOLibkky/
addpath ../utils/

warning off;

% choose function
numExperiments = 2;
funcChoice = 4;

% Function 1
if funcChoice == 1
  numDims = 5;
  func = @(t) t.^2 /4 + sin(7*t);
  bounds = [0 1];
  maxNumFunEvals = 15;

% Function 2
elseif funcChoice == 2
  numDims = 20;
  [func, funcProperties] = get3ModalFunction(numDims);
  bounds = funcProperties.bounds;
  maxNumFunEvals = 100; min(8 * 2^numDims, 1000);
  maxPt = funcProperties.maxPt;
  maxVal = func(funcProperties.maxPt);
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

elseif funcChoice == 4
  numDims = 20;
  numDimsPerGroup = 4;
  [func, funcProperties] = getAdditiveFunction(numDims, numDimsPerGroup);
  bounds = funcProperties.bounds;
  maxNumFunEvals = min(8 * 2^numDims, 1000);
  maxPt = funcProperties.maxPt;
  maxVal = funcProperties.maxVal;
  fprintf('True maxVal, maxPt: %f, %s\n', maxVal, mat2str(maxPt) );
end

% Set parameters for Bayesian Optimization
boParams.optPtStdThreshold = 0.002;
boParams.alBWLB = 1e-5;
boParams.alBWUB = 5;
boParams.numInitPts = 4*numDims;
boParams.gpNoiseLevel = 0.5; % * std(f(th));

% Initialize arrays for storing the history
totalNumQueries = maxNumFunEvals + boParams.numInitPts;
boHistories = zeros(numExperiments, totalNumQueries);
ucbHistories = zeros(numExperiments, totalNumQueries);
rpBoHistories = zeros(numExperiments, totalNumQueries);
rbBoHistories = zeros(numExperiments, totalNumQueries);
randHistories = zeros(numExperiments, totalNumQueries);

for expIter = 1:numExperiments

  % 3.1 Call BO
  boParams.utilityFunc = 'EI';
  [boMaxVal, boMaxPt, boQueries, boVals, boHist] = ...
    bayesOptGP(func, bounds, maxNumFunEvals, boParams);
  boHistories(expIter, :) = boHist';

  % Call BO with UCB
  boParams.utilityFunc = 'UCB';
  [ucbMaxVal, ucbMaxPt, ucbQueries, ucbVals, ucbHist] = ...
    bayesOptGP(func, bounds, maxNumFunEvals, boParams);
  ucbHistories(expIter, :) = ucbHist';

%   % 3.2 Call RP Bayes
%   projDims = 1;
%   [rpMaxVal, rpMaxPt, rpBoQueries, rpBoVals, rpBoHist] = ...
%     boRPGP(func, bounds, projDims, maxNumFunEvals, boParams);
%   rpBoHistories(expIter, :) = rpBoHist;

  % 3.3 Call RandomBases
%   [rbMaxVal, rbMaxPt, rbBoQueries, rbBoVals, rbBoHist] = ...
%     boRandomBasesGP(func, bounds, maxNumFunEvals, boParams);
%   rbBoHistories(expIter, :) = rbBoHist;

  % 3.4 Do Random calls
  randQueries = bsxfun(@plus, ...
          bsxfun(@times, rand(totalNumQueries, numDims), ...
                 (bounds(:,2) - bounds(:,1))' ), bounds(:,1)');
  randVals = func(randQueries);
  randHistories(expIter, :) = max(randVals(cumsum(triu(ones(length(randVals))))));

end

  % finally Call DIRECT
  diRectOpts.maxevals = totalNumQueries;
  [~, ~, diRectHist] = diRectWrap(func, bounds, diRectOpts);
  numFuncEvaluations = diRectHist(end, 2);

  boHistMean = mean(boHistories, 1); boHistStd = std(boHistories, 1);
  ucbHistMean = mean(ucbHistories, 1); ucbHistStd = std(ucbHistories, 1);
  rpBoHistMean = mean(rpBoHistories, 1); rpBoHistStd = std(rpBoHistories, 1);
  rbBoHistMean = mean(rbBoHistories, 1); rbBoHistStd = std(rbBoHistories, 1);
  randHistMean = mean(randHistories, 1); randHistStd = std(randHistories, 1);

  % 4. Plot iteration statistics
  figure;  hold on,
  loglog(1:totalNumQueries, boHistMean, 'b'); hold on,
  loglog(1:totalNumQueries, ucbHistMean, 'c'); hold on,
  loglog(1:totalNumQueries, rpBoHistMean, 'g'); 
  loglog(1:totalNumQueries, rbBoHistMean, 'm'); 
  loglog(diRectHist(:,2),diRectHist(:,3), 'r');
  loglog(1:totalNumQueries, randHistMean, 'k');
  legend('EI', 'UCB', 'randomProjectBO', 'randomBaseBO', 'DiRect', 'Rand');
  if numExperiments > 1
  errorbar(1:totalNumQueries, boHistMean, boHistStd/sqrt(numExperiments), 'Color', 'b');
  errorbar(1:totalNumQueries, ucbHistMean, ucbHistStd/sqrt(numExperiments), 'Color', 'c');
  errorbar(1:totalNumQueries, rpBoHistMean, rpBoHistStd/sqrt(numExperiments), 'Color', 'g');
  errorbar(1:totalNumQueries, rbBoHistMean, rbBoHistStd/sqrt(numExperiments), 'Color', 'm');
  errorbar(1:totalNumQueries, randHistMean, randHistStd/sqrt(numExperiments), 'Color', 'k');
  end
  loglog([0, totalNumQueries], [maxVal, maxVal], 'k')';
%   axis([20*numDims totalNumQueries*1.1, maxVal - (maxVal-randHistMean(1))/1, maxVal*1.1]);

  % 5. Plot the function
  if numDims == 1
    figure;
    res = 200;
    th = linspace(0,1, res);
    % Plot
    plot(th, func(th), 'b'); hold on,
    plot(boQueries, boVals, 'k.');
    plot(maxPt, maxVal, 'rx', 'MarkerSize', 10);
    plot(rpMaxPt, rpMaxVal, 'g*', 'MarkerSize', 10);
    plot(rbMaxPt, rbMaxVal, 'm*', 'MarkerSize', 10);
  elseif numDims == 2
    figure;
    plot2DFunction(@(t) (func(t)), bounds, 'contour'); hold on,
    plot(maxPt(1), maxPt(2), 'rx', 'MarkerSize', 10, 'LineWidth', 3);
    plot(rpMaxPt(1), rpMaxPt(2), 'g*', 'MarkerSize', 10, 'LineWidth', 3);
    plot(rbMaxPt(1), rbMaxPt(2), 'm*', 'MarkerSize', 10, 'LineWidth', 3);
    title('func');
  end

