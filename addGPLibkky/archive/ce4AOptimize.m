function [A, optVal] = ce4AOptimize(func, d, p, potInitAs, numInitPts)
% Okay, I am clearly running out of ideas for function names now ..
% This method combines the gradient method of Wen, Yin along with a cross
% entropy method to avoid local maximum.

  % Prelims
  if isempty(numInitPts), numInitPts = d*p;
  end
  if isempty(potInitAs), numPotentialInitAs = 0;
  else numPotentialInitAs = size(potInitAs, 3);
  end
  numElitePts = 5;
  numInitGDIters = 100; 
  numFinalGDIters = 300;
%   numInitGDIters = 5; 
%   numFinalGDIters = 5; 30;

  % Set optimization parameters
  stiefelOpts.record = 0;
  stiefelOpts.xtol = 1e-5;
  stiefelOpts.gtol = 1e-5;
  stiefelOpts.ftol = 1e-5;

  % Initialization
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % First initialize
  initAs = zeros(d, p, numInitPts);
  initAs(:,:,1:numPotentialInitAs) = potInitAs; 
  for i = (numPotentialInitAs+1):numInitPts
    [Q,~] = qr(rand(d,p), 0);
    initAs(:,:,i) = Q;
  end

  initVals = zeros(numInitPts, 1);
  parfor initIter = 1:numInitPts
    initVals(initIter) = func( initAs(:,:,initIter) );
  end

  % TODO from here
  %%%%%%%%%%%%%%%%%%%%%

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

