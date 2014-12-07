function [A, optVal] = decompOptimize(func, D, p, potInitAs, numInitPts)
% Optimizes to find the best decomposition


  % Prelims
  if isempty(numInitPts), numInitPts = 15*D*p;
  end
  if isempty(potInitAs), numPotentialInitAs = 0;
  else numPotentialInitAs = size(potInitAs, 3);
  end
  numElitePts = 1;
  numGDIters = 100;
  numInitGDIters = 4; 
  numDecompIters = 1;

  % Set optimization parameters
  stiefelOpts.record = 0;
  stiefelOpts.xtol = 1e-5;
  stiefelOpts.gtol = 1e-5;
  stiefelOpts.ftol = 1e-5;

  % Initialization
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % First initialize
  initAs = zeros(D, p, numInitPts);
  initAs(:,:,1:numPotentialInitAs) = potInitAs; 
  for i = (numPotentialInitAs+1):numInitPts
    [Q,~] = qr(rand(D,p), 0);
    initAs(:,:,i) = Q;
  end

  % Now perform numInitGDIters iterations on them and project to permutation
  % matrices. Then compute the objective.
  initVals = zeros(numInitPts, 1);
  stiefelOpts.mxitr = numInitGDIters;
  for i = 1:numInitPts
    Ai = initAs(:,:,i);
    Ai = OptStiefelGBB(Ai, func, stiefelOpts);
    Ai = orthToPermutation(Ai);
    initAs(:,:,i) = Ai;
    initVals(i) = func( initAs(:,:,i) );
  end

  % Now pick the best numElitePts matrices
  [~, uniqueIdxs] = unique(initVals);
  if numel(uniqueIdxs) < numElitePts, eliteIdxs = uniqueIdxs;
  else eliteIdxs = uniqueIdxs(1:numElitePts);
  end
  eliteAs = initAs(:,:,eliteIdxs);
  numElitePts = numel(eliteIdxs);
  
  % now, do this iteratively.
  eliteVals = zeros(numElitePts, 1);
  stiefelOpts.mxitr = numGDIters;
  for j = 1:numDecompIters
    for i = 1:numElitePts
      Ai = eliteAs(:,:,i);
      Ai = OptStiefelGBB(Ai, func, stiefelOpts);
      [Ai, permi] = orthToPermutation(Ai);
%       fprintf('i,j= %d,%d', i,j); Ai,
      eliteAs(:,:,i) = Ai;
      eliteVals(i) = func(Ai);
%       fprintf('i,j,eliteVals(i)= %d,%d,%f, %s\n', i,j,eliteVals(i), mat2str(permi));
    end
  end

  % Pick the best
  [optVal, bestPt] = min(eliteVals);
  A = eliteAs(:,:,bestPt);

end

