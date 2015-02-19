% This performs a few tests for kde.m

clear all;
close all;
N = 300;

% Test 1: Gaussian 1D example
% ===========================

fprintf('Test 1: Gaussian 1D\n');
X = randn(N, 1);
[~, f, bandwidth] = kde(X);
% plot the curve and estimate the density. Verify that it integrates to 1.
x = linspace(-6,6,1000)';
p = f(x);
plot(x, p, 'b', x, normpdf(x, 0, 1), 'r'); hold on,
plot(X, 0.1*rand(size(X)), 'kx');
hold off;
titlestr = sprintf('Estimated(b) and True(r)\nh = %f', bandwidth);
title(titlestr);
area = numerical_1D_integration(x, p);
fprintf('Bandwidth : %f, Area under estimated curve: %f\n', bandwidth, area);
fprintf('Press any key to continue ...\n\n');
pause;


% Test 2: 2D Gaussian Mixture Model
% =================================
  
fprintf('Test 2: GMM 2D\n');
p1 = 0.7;
p2 = 1 - p1;
mu1 = [0; 0];
S1 = [1 0.5; 0.5 2];
mu2 = [1; 3];
S2 = [3 -1; -1 4];
Z = double(rand(N,1) < p1);
X = [Z, Z] .* bsxfun(@plus, randn(N, 2) * chol(S1), mu1') + ...
    [1-Z, 1-Z] .* bsxfun(@plus, randn(N, 2) * chol(S2), mu2');
fprintf('performing KDE\n');
[~, f, bandwidth] = kde(X);
fprintf('done KDE\n');
% Pl0t the density.
t = linspace(-10, 10, 1000)'; [T1, T2] = meshgrid(t, t); T =[T1(:), T2(:)];
p = f(T); P = reshape(p, 1000, 1000);
true_density = p1 * mvnpdf(T, mu1', S1) + p2 * mvnpdf(T, mu2', S2);
True_Density = reshape(true_density, 1000, 1000);
contour(T1, T2, P, 'color', 'b'); hold on,
contour(T1, T2, True_Density, 'color', 'r');
plot(X(:,1), X(:,2), 'kx');
titlestr = sprintf('Estimated(b) vs True(r)\nh = %f', bandwidth);
title(titlestr);
area = numerical_2D_integration(P, T1, T2);
fprintf('Bandwidth : %f, Area under estimated curve: %f\n', bandwidth, area);

