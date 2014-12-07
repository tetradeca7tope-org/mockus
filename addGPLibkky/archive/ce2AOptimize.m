function [A, optVal] = ce2AOptimize(func, d, p, fineOptimize)
% Okay, I am clearly running out of ideas for function names now ..
% This method combines the gradient method of Wen, Yin along with a cross
% entropy method to avoid local maximum.

  if ~exist('fineOptimize', 'var')
  % fineOptimize tells whether to do some exhaustive optimization at the end.
    fineOptimize = true;
  end

  % Prelims
  numSpawnPts = 1000;
  numElitePts = 5;
  numNonElitePts = 0; % Just to avoid local optima
  numAnchorPts = numElitePts + numNonElitePts;
  numSamples = numAnchorPts * numSpawnPts;
  numGDItersPerStep = 2; 
  numCEIters = 25;
  lambda = 0.2/(d*p); 

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
    % First go through each sample and compute the value
    parfor sampleIter = 1:numSamples
      Ai = currAs(:,:,sampleIter);
      sampleVals(sampleIter) = func(Ai);
    end

    % Now Pick the new set of anchor points 
    [~,sortedIdxs] = sort(sampleVals);
    eliteSamples = sortedIdxs(1:numElitePts);
    remSamples = sortedIdxs(numElitePts+1:end);
    nonEliteSamples = remSamples( randsample(numSamples-numElitePts, ...
                                            numNonElitePts));
    anchorSamples = [eliteSamples; nonEliteSamples];
    anchorAs = currAs(:,:,anchorSamples);
    % Save the best guy
    currBestSample = sortedIdxs(1);
    sampleVals(anchorSamples)', currBestSample,
    if optVal > sampleVals(currBestSample)
      optVal = sampleVals(currBestSample);
      A = currAs(:,:,currBestSample);
    end

    if ceIter ~= numCEIters
      % Perform Gradient Descent on the anchor points
      for anchorIter = 1:numAnchorPts
        Ai = anchorAs(:,:,anchorIter); 
        [Ai, out] = OptStiefelGBB(Ai, func, stiefelOpts);
        anchorAs(:,:,anchorIter) = Ai;
      end

      % Now create the new set of samples
      newAs = zeros(d,p,numSamples);
      for i = 1:numAnchorPts
        newAs(:,:, (i-1)*numSpawnPts + 1) = anchorAs(:,:,i);
        for j = 1:(numSpawnPts-1)
          B = (1-lambda) * currAs( anchorSamples(i) ) + lambda*rand(d,p);
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
    for i = 1:10*numCEIters
      [An, out] = OptStiefelGBB(An, func, stiefelOpts);
      outVal = func(An);
      if optVal > outVal 
        optVal = outVal;
        A = An; 
      end
    end
  end
  
end

