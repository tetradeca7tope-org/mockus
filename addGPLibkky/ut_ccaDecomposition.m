% Unit Test for CCA Decomposition

clear all;
close all;
clc; clc;

% First get an additive function
numDims = 8;
numDimsPerGroup = 1;
% [f, fProps] = getRotAdditiveFunction(numDims, numDimsPerGroup);
[f, fProps] = getAdditiveFunction(numDims, numDimsPerGroup); fProps.A = eye(numDims);

numMCMCSamples = 10000 * numDims;
mcmcPropStd = 0.1 * sqrt(numDims);

% Now do CCA
[A, queries, vals, samples] = ccaDecomposition(f, numDims, numMCMCSamples, ...
  mcmcPropStd);

% Print out As
fProps.A, A, 
%inv(A)'*fProps.A,
err = norm(abs(fProps.A) - abs(A)),
err = norm(abs(fProps.A) - abs(A')),
fprintf('Error: %f\n', err);

plot(samples(:,1), samples(:,2), 'rx');

if numDims == 2
  figure;
  bounds = [-1 1; -1 1];
  plot2DFunction(f, bounds, 'contour');
end

