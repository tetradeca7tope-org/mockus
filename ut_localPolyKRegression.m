% Unit test for locally polynomial kernel regression

close all;
clear all;

% 1 dimensional test
num_train_pts = 10;
resolution = 100;
h = 0.05;

% prep the function handle
f = @(arg) arg.^2 + 0.3*arg + 0.2 - sin(6*pi*arg); % true function

% prep the data
X = linspace(0, 1, num_train_pts/2 + 2)'; X = X(2:end-1);
X = [X; rand(num_train_pts/2, 1)];
Y = f(X);
th = linspace(0, 1, resolution)';

% TEST for localPolyKRegression.m
% ===============================
fprintf('1D test for localPolyKRegression \n');
% fit the function for different poly orders
est_th0 = localPolyKRegression(th, X, Y, h, 0);
est_th1 = localPolyKRegression(th, X, Y, h, 1);
est_th2 = localPolyKRegression(th, X, Y, h, 2);
% est_th5 = locallyPolyKernelRegression(th, X, Y, h, 5);
plot(X, Y, 'kx'); hold on
plot(th, est_th0, 'r-x');
plot(th, est_th1, 'r-.');
plot(th, est_th2, 'g--');
% plot(th, est_th5, 'r-o');
plot(th, f(th), 'b');
legend('pts', 'k = 0', 'k = 1', 'k = 2', 'k=5', 'true');
title('Test for localPolyKRegression');

% TEST for localPolyKRegressionCV.m
% =================================
fprintf('\n1D test for localPolyKRegressionCV \n');
[est_th, best_h, best_k] = localPolyKRegressionCV(th, X, Y);
figure;
plot(X, Y, 'kx'); hold on
plot(th, f(th), 'b');
plot(th, est_th, 'r-');
legend('pts', 'true', 'estimated');
title('Test for localPolyKRegressionCV');
axis([0 1 -1 2.5]);
fprintf('Chosen (h, k) = (%0.4f, %d)\n', best_h, best_k);

