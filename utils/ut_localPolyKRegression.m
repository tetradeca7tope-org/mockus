% Unit test for locally polynomial kernel regression

close all;
clear all;

% 1 dimensional test
num_train_pts = 1000;
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
tic,
plot(th, est_th0, 'r-x');
plot(th, est_th1, 'r-.');
plot(th, est_th2, 'g--');
toc,
% plot(th, est_th5, 'r-o');
plot(th, f(th), 'b');
legend('pts', 'k = 0', 'k = 1', 'k = 2', 'true');
title('Test for localPolyKRegression');

% TEST for localPolyKRegressionCV.m
% =================================
fprintf('\n1D test for localPolyKRegressionCV \n');
tic,
[est_th, best_h, best_k] = localPolyKRegressionCV(th, X, Y);
toc,
figure;
plot(X, Y, 'kx'); hold on
plot(th, f(th), 'b');
plot(th, est_th, 'r-');
legend('pts', 'true', 'estimated');
title('Test for localPolyKRegressionCV');
axis([0 1 -1 2.5]);
fprintf('Chosen (h, k) = (%0.4f, %d)\n', best_h, best_k);

% TEST for localPolyKRegressionCV.m in 2D
% =======================================
fprintf('\n2D test for localPolyKRegressionCV \n');

% Define the function and plot it.
figure;
f = @(X) (X*[1; 1.5] - 1.25).^2 + sin(X*[2; 1.3]);
t = linspace(0, 1, resolution); [T1, T2] = meshgrid(t, t);
th = [T1(:), T2(:)];
thf = f(th);
fT = reshape(thf, resolution, resolution);
mesh(T1, T2, fT); hold on,

% Define the inputs
X = rand(num_train_pts, 2);
Y = f(X);

% perform the regression
tic,
[est_th, best_h, best_k] = localPolyKRegressionCV(th, X, Y);
toc,
EstTh = reshape(est_th, resolution, resolution);
mesh(T1, T2, EstTh);

plot3(X(:,1), X(:,2), Y, 'kx');
title('2D Test for localPolyKRegressionCV');
fprintf('Chosen (h, k) = (%0.4f, %d)\n', best_h, best_k);

