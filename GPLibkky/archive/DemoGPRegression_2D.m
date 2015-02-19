% A Demo for GP regression.
close all;
clear all;

% runtime params
NOISE_LEVEL = 0.15;

% function
f = @(X1, X2) (sin(4*X1 + 2*X2) + (X1 + 1.4*X2).^2 + X2);

% Generate data
x1 = [0.2, 0.4, 0.6, 0.8]; x2 = [0.3 0.5 0.6 0.9];
[X1,X2] = meshgrid(x1,x2); X = [X1(:), X2(:)];
Y = f(X(:,1), X(:,2));

sigmaSm = 0.5
sigmaPr = 2*std(Y);

% Generate test points
N = 100;
x1 = linspace(0,1,N)'; x2 = linspace(0,1,N)';
[X1,X2] = meshgrid(x1,x2); Xte = [X1(:), X2(:)];
Ytrue = f(X1, X2);

% Perform GP-Regression
hyperParams.sigmaSm = sigmaSm;
hyperParams.sigmaPr = sigmaPr;
hyperParams.noise = NOISE_LEVEL*ones(size(X,1),1);;
hyperParams.meanFunc = [];
runtimeParams.plotOn = false;
[post_mean, ~, K] = GPRegression(X, Y, Xte, hyperParams, runtimeParams);

% Plot the results
PLOTFUNC = @surface;
Yte = reshape(post_mean, N, N);
figure; PLOTFUNC(X1, X2, Yte ); hold on
plot3(X(:,1), X(:,2), Y, 'rx', 'LineWidth', 4, 'MarkerSize', 10);
figure; PLOTFUNC(X1, X2, Ytrue);
fprintf('Avg sum of squared errors: %f\n', mean(mean( (Yte - Ytrue).^2 )));

figure;
% Now draw some samples
num_samples = 100;
gp_samples = GPDrawSamples(post_mean, K, num_samples);
% plot the samples
for i = 1:num_samples
  curr_sample = reshape(gp_samples(i,:), N, N);
  PLOTFUNC(X1, X2, Yte ); hold on,
%   PLOTFUNC(X1, X2, curr_sample, 'Color', 'm'); hold off,
  PLOTFUNC(X1, X2, curr_sample, 'FaceLighting', 'gouraud'); hold off,
%   % Also print out the avg likelihood of this function
%   fprintf('sample %d: avg-log-likelihood = %f\n', ...
%           i, GPAvgLogLikelihood(post_mean, K, gp_samples(i,:)'));
  pause,
end

