function [A, optVal] = ceAOptimize(func, d, p, fineOptimize)
% Okay, I am clearly running out of ideas for function names now ..
% This method combines the gradient method of Wen, Yin along with a cross
% entropy method to avoid local maximum.

  if ~exist('fineOptimize', 'var')
  % fineOptimize tells whether to do some exhaustive optimization at the end.
    fineOptimize = true;
  end

  % Prelims
  numSpawnPts = 10;
  numElitePts = 10;
  numNonElitePts = 0; % Just to avoid local optima
  numAnchorPts = numElitePts + numNonElitePts;
  numSamples = numAnchorPts * numSpawnPts;
  numGDItersPerStep = 2; 
  numCEIters = 25;
  lambda = 0.9; % perturbation is 1 - lambda

  % Set optimization parameters
  stiefelOpts.record = 0;
  stiefelOpts.mxitr = numGDItersPerStep;
  stiefelOpts.xtol = 1e-5;
  stiefelOpts.gtol = 1e-5;
  stiefelOpts.ftol = 1e-5;


  % First initialize
  currAs = zeros(d, p, numSamples);
  for i = 1:numSamples
    currAs(:,:,i) =  orth(randn(d,p));
  end
  
  optVal = inf;
  for ceIter = 1:numCEIters
    
    sampleVals = zeros(numSamples, 1);
    % First go through each sample and ascend in a few directions.
    for sampleIter = 1:numSamples
      Ai = currAs(:,:,sampleIter);
%       Ai, Ai'*Ai, 
      [Ai, out] = OptStiefelGBB(Ai, func, stiefelOpts);
      sampleVals(sampleIter) = func(Ai);
      currAs(:,:,sampleIter) = Ai;
    end

    % Now Pick the new set of non-Elite samples
    [~,sortedIdxs] = sort(sampleVals);
    eliteSamples = sortedIdxs(1:numElitePts);
    remSamples = sortedIdxs(numElitePts+1:end);
    nonEliteSamples = remSamples( randsample(numSamples-numElitePts, ...
                                            numNonElitePts));
    anchorSamples = [eliteSamples; nonEliteSamples];

    % Save the best guy
    sampleVals(anchorSamples)',
    currBestSample = sortedIdxs(1);
    currBestSample,
    if optVal > sampleVals(currBestSample)
      optVal = sampleVals(currBestSample);
      A = currAs(:,:,currBestSample);
    end

    if ceIter ~= numCEIters
      % Now create the new set of samples
      newAs = zeros(d,p,numSamples);
      for i = 1:numAnchorPts
        newAs(:,:, (i-1)*numSpawnPts + 1) = currAs(:,:, anchorSamples(i) );
        for j = 1:(numSpawnPts-1)
          B = lambda * currAs( anchorSamples(i) ) + (1 - lambda)*randn(d,p);
          [B,~] = qr(B,0);
          newAs(:,:, (i-1)*numSpawnPts + j + 1) = B;
%           (i-1)*numSpawnPts + j + 1, newAs(:,:, (i-1)*numSpawnPts + j + 1), pause,
        end
      end
      % Finally save it back
      currAs = newAs;
    end
    
    fprintf('ceIter: %d, optVal: %0.4f\n', ceIter, optVal);
    A,
  end

  % Finally do some descent on the best point again.
  if fineOptimize
    An = A;
    for i = 1:numCEIters
      [An, out] = OptStiefelGBB(An, func, stiefelOpts);
      outVal = func(An);
      if optVal > outVal 
        optVal = outVal;
        A = An; 
      end
    end
  end
  
end

