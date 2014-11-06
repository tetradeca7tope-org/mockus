% Unit test for genPolyFeatures and genDecompPolyFeatures

clc;
X = [1:4];
X = [X; 2*X];

fprintf('genPolyFeatures.m\n');
[F, G] = genPolyFeatures(X, 1),

fprintf('genDecompPolyFeatures.m\n');
d = 5;
X = [1:d; 2:1:(d+1)];
A = eye(d); [F, G] = genDecompPolyFeatures(X, A, 2, 1), full(F), full(G),
A = randn(5); A = orth(A); [F, G] = genDecompPolyFeatures(X, A, 2, 1), full(F), full(G),

size(F),
size(G),
