% Unit test for spectralCluster.m
% Uses the 2-moons dataset.

% Params
h = 0.05; % Bandwidth for Gaussian Kernel

% Generate 2 moon dataset
N = 100;
theta = linspace(0, pi, N)';
X1 = [0.37 + 0.27*cos(theta), 0.45 + 0.27*sin(theta)];
X1 = X1(randperm(N), :);
X2 = [0.64 + 0.27*cos(theta), 0.55 - 0.27*sin(theta)];
X2 = X2(randperm(N), :);
X = [X1; X2];
% Add some noise.
X = X + 0.01*randn(size(X));
% X = X(randperm(size(X,1)), :);

% Construct similarity matrix
D = dist2(X, X);
A = exp(-D/(2*h^2) );

% Perform Spectral Clustering
labels = spectralCluster(A, 2);
Z1 = X(labels ==1, :);
Z2 = X(labels ==2, :);
figure;
plot(Z1(:,1), Z1(:,2), 'bx', Z2(:,1), Z2(:,2), 'ro');

