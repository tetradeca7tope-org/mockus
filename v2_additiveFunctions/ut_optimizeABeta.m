% Unit Test for optimizeABeta.m

n = 100;
D = 10;
% f = @(X) (X(:,1) + 2*X(:,2)).^2 + 0.5 * (X(:,3) + 0.8* X(:,4)).^3 + X(:,4).^4;
% f = @(X) (X(:,1) + 2*X(:,2) + X(:,3) + 0.8* X(:,4)).^3 + X(:,4).^4;(X(:,1) +
f = @(X) (X(:,1) + X(:,2)).^2 + ( X(:,3) + 2*X(:,4)).^2 + X(:,3).^3;


X = rand(100, D);
Y = f(X);

params.verbose = true;
params.maxOptAIters = 10;
params.numInits = 10;
params.numIters = 10;
params.d = 2;
params.lambda = 0.1;
[beta, A] = optimizeABeta(Y, X, params);
