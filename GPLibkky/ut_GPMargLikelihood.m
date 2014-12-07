% A Demo for GP regression with Leave One Out Cross Validation.
close all;
clear all;
addpath ../utils/
addpath ../ancillary/

% Test 1
% ==============================================================================
fprintf('Test 1\n======================\n');
% Generate data
m = 5; X = linspace(0,1,m)';
f = @(X) X.^2 + X.*sin(4*X) + 0*randn(size(X,1),1);
y = f(X);
% X = [-4 -3 -1 0 1.5]'; X = (X + 4)/5.5; y = [-2 0 1 2 -1]'; 

% Generate test points
N = 100;
Z = linspace(0,1,N)';

% Perform GP-Regression
% hyperparams.noise = 0.1 * ones(size(X,1),1);
hyperparams.meanFunc = [];
hyperparams.sigmaSmRange = [];
hyperparams.sigmaPrRange = [0.1 10] * std(y);
% hyperparams.sigmaSm = 0.24793;
% hyperparams.sigmaPr = 2.99450;
hyperparams.sigmaSm = 0;
hyperparams.sigmaPr = 0;
[post_mean, K, funcH, sigmaSm, sigmaPr] = ...
  GPMargLikelihood(X, y, Z, hyperparams);

  % Now draw some samples
  num_samples = 100;
  gp_samples = GPDrawSamples(post_mean, K, num_samples);
  % plot the samples
  figure;
  hold on,
  for i = 1:num_samples
    plot(Z, gp_samples(i,:), 'm-');
  end
plot(Z, post_mean, 'b', 'LineWidth', 2); hold on,
plot(Z, post_mean + diag(K), 'b--', 'LineWidth', 1); hold on,
plot(Z, post_mean - diag(K), 'b--', 'LineWidth', 1); hold on,
plot(Z, f(Z), 'g', 'LineWidth', 2);
plot(X, y, 'kx', 'MarkerSize', 10, 'Linewidth', 2);
fprintf('Chosen: sigmaSm: %f, sigmaPr: %f, \nlikl of post-mean: %f\n', ...
  sigmaSm, sigmaPr, GPAvgLogLikelihood(post_mean, K, post_mean));
% pause;
% fprintf('Paused .... \n');

% Test 2
% ==============================================================================
fprintf('Test 2\n======================\n');

% function
f = @(X1, X2) ( sin(4*X1 + 2*X2) + (X1 + 1.4*X2).^2 + X2);

% Generate Data
x1 = [0.2 0.4 0.6 0.8]; x2 = [0.3 0.5 0.6 0.9];
x1 = linspace(0.1, 0.9, 5); x2 =x1;
[X1, X2] = meshgrid(x1, x2); X = [X1(:), X2(:)];
Y = f(X(:,1), X(:,2));

% Generate test points
N = 50;
x1 = linspace(0,1,N)';
[X1, X2] = meshgrid(x1, x1); Xte = [X1(:), X2(:)];
Ytrue = f(X1, X2);

% Perform GP-Regression
% hyperparams.noise = NOISE_LEVEL*ones(size(X,1),1);
hyperparams.meanFunc = [];
hyperparams.sigmaSmRange = [];
hyperparams.sigmaPrRange = [];
hyperparams.sigmaSm = 0;
hyperparams.sigmaPr = 0;
[post_mean, K, funcH, sigmaSm, sigmaPr] = ...
  GPMargLikelihood(X, Y, Xte, hyperparams);

% Plot the results
figure;
PLOTFUNC = @surface;
Yte = reshape(post_mean, N, N);
PLOTFUNC(X1, X2, Yte ); hold on
plot3(X(:,1), X(:,2), Y, 'rx', 'LineWidth', 4, 'MarkerSize', 10);
title('estimate');

figure; 
plot3(X(:,1), X(:,2), Y, 'rx', 'LineWidth', 4, 'MarkerSize', 10);
PLOTFUNC(X1, X2, Ytrue);
title('truth');
fprintf('Avg sum of squared errors: %f\n', mean(mean( (Yte - Ytrue).^2 )));


% Test 3 : Higher Dimensions
% ==============================================================================
fprintf('Test 3\n======================\n');

numDims = 40;
numTrainData = numDims^2;
numTestData = numDims^2;
f = get2ModalFunction(numDims);
XTrain = rand(numTrainData, numDims);
YTrain = f(XTrain);
XTest = rand(numTestData, numDims);
YTest = f(XTest);

hyperparams.meanFunc = [];
hyperparams.sigmaSmRange = [];
hyperparams.sigmaPrRange = [];
hyperparams.sigmaSm = 0;
hyperparams.sigmaPr = 0;
[postMean] = GPMargLikelihood(XTrain, YTrain, XTest, hyperparams);

% Report error
testError = norm(postMean - YTest)/sqrt(numTestData);
fprintf('Error: %0.5f\n', testError);
fprintf('Range: %0.5f\n', max([YTrain; YTest]) - min([YTrain; YTest]));
fprintf('Std: %0.5f\n', std([YTrain; YTest]));

