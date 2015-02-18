% A Demo for GP regression with Leave One Out Cross Validation.
close all;
clear all;

% Generate data
m = 5000; X = linspace(0,1,m)'; y = X.^2 + X.*sin(4*X) + 0*randn(m,1);
% X = [-4 -3 -1 0 1.5]'; X = (X + 4)/5.5; y = [-2 0 1 2 -1]'; 

candidates.sigmaSmVals = logspace(-0.5, 1, 8)' * std(X);
candidates.sigmaPrVals = logspace(-0.5, 1, 6)' * std(y);

% Generate test points
N = 100;
Z = linspace(0,1,N)';

% Perform GP-Regression
hyperparams.noise = zeros(size(X,1),1);
hyperparams.meanFunc = [];
[post_mean, K, sigmaSm, sigmaPr] = ...
  GPKFoldCV(X, y, Z, 10,candidates, hyperparams);

% Now draw some samples
% num_samples = 100;
% gp_samples = GPDrawSamples(post_mean, K, num_samples);
% % plot the samples
% figure;
% hold on,
% for i = 1:num_samples
%   plot(Z, gp_samples(i,:), 'm-');
% end
% plot the mean on top of it all
% pause,
plot(Z, post_mean, 'b', 'LineWidth', 2); hold on,
plot(X, y, 'kx');
fprintf('Chosen: sigmaSm: %f, sigmaPr: %f, \nlikl of post-mean: %f\n', ...
  sigmaSm, sigmaPr, GPAvgLogLikelihood(post_mean, K, post_mean));

