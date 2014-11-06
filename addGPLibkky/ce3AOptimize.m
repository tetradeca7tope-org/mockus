function [A, optVal] = ce2AOptimize(func, d, p, fineOptimize)
% Okay, I am clearly running out of ideas for function names now ..
% This method combines the gradient method of Wen, Yin along with a cross
% entropy method to avoid local maximum.

  if ~exist('fineOptimize', 'var')
  % fineOptimize tells whether to do some exhaustive optimization at the end.
    fineOptimize = true;
  end

  % Prelims
  numInitPts = 10*d*p;
  numElitePts = 5;
  numInitGDIters = 100; 
  numFinalGDIters = 300;
  numGDIters = 200;

  % Set optimization parameters
  stiefelOpts.record = 0;
  stiefelOpts.xtol = 1e-5;
  stiefelOpts.gtol = 1e-5;
  stiefelOpts.ftol = 1e-5;

  % First initialize
  initAs = zeros(d, p, numInitPts);
  for i = 1:numInitPts
    [Q,~] = qr(rand(d,p), 0);
%     [Q,~] = qr(randn(d,p), 0);
    initAs(:,:,i) = Q;
  end
%   initAs(:,:,1) = eye(d,p);

  initVals = zeros(numInitPts, 1);
  parfor initIter = 1:numInitPts
    initVals(initIter) = func( initAs(:,:,initIter) );
  end

  % Now Pick the new set of anchor points 
  [~, sortedIdxs] = sort(initVals);
  elitePts = sortedIdxs(1:numElitePts);
  eliteAs = initAs(:,:,elitePts);
  eliteVals = initVals(elitePts);
%   eliteVals,

  stiefelOpts.mxitr = numInitGDIters;
  for i = 1:numElitePts
  % Maximize each of these points for a few iterations
    Ai = eliteAs(:,:,i);
    Ai = OptStiefelGBB(Ai, func, stiefelOpts);
    eliteAs(:,:,i) = Ai;
    eliteVals(i) = func(Ai);
  end
%   eliteVals,

  % Pick the best
  [~, bestPt] = min(eliteVals);
  bestA = eliteAs(:,:,bestPt);

  % Now maximize over the best
  stiefelOpts.mxitr = numFinalGDIters;
  A = OptStiefelGBB(bestA, func, stiefelOpts);
  optVal = func(A);
  
end

