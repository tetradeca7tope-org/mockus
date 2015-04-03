function [maxVal, maxPt, boQs, boVals, chosen, history] = ...
bayesOptDecideAddGP(oracle, decomp, bounds, numIters, params, numDims)

% decomp is a cell array, each cell is a struct with the following 
% fields: d, M. 

DECOMP_KNOWN = 'known';
DECOMP_LEARN = 'learn';
DECOMP_RAND = 'random';
DECOMP_PLEARN = 'partialLearn';
DECOMP_DECIDE = 'decide';

if ~strcmp(params.decompStrategy, DECOMP_DECIDE)
  [maxVal, maxPt, boQueries, boVals, history] = ...
  bayesOptDecompAddGP(oracle, decomp, bounds, numIters, params);

else
  numDecompChoices = numel(decomp);
  numOutIters = numDecompChoices;
  numInnerIters = floor(numIters / numOutIters);
  
  decompNext.d = decomp{1}.d;
  decompNext.M = decomp{1}.M;
  
  boQs = [];
  boVals = [];
  
  % Initialize array that saves chosen (d,M) at each stage
  chosen = nan(numOutIters,2);
  
  for cout = 1:numOutIters
    chosen(cout,:) = [decompNext.d; decompNext.M];
    
    params.decompStrategy = 'decide';
    [maxVal, maxPt, thisBoQs, thisBoVals, thisHistory, gpHyperParams]=...
      bayesOptDecompAddGP(oracle, decompNext, bounds, numInnerIters, params);
    
    if (cout == 1)
      boQs = [boQs; thisBoQs];
      boVals = [boVals; thisBoVals];
    else
      len = numel(thisBoVals);
      boQs = [boQs; thisBoQs(len-numInnerIters+1:end,:)];
      boVals = [boVals; thisBoVals(len-numInnerIters+1:end)]; 
    end

    dummyPts = zeros(0,numDims);

    currBestMll = inf;
    currBestDecompIdx = -1;

    for i=1:numDecompChoices
      testDecomp = decomp{i};
      numGroups = testDecomp.M;

      % gpHyperParams.sigmaSmRanges
      
       [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, learnedDecomp, mll] = ...
       addGPDecompMargLikelihood(boQs, boVals, dummyPts, testDecomp, gpHyperParams);
      
    
      if mll < currBestMll 
        currBestMll = mll;
        currBestDecompIdx = i;
      end
    end
    
    decompNext = decomp{currBestDecompIdx}; 
    fprintf('\nChoose (d,M) = (%d,%d)\n', decompNext.d, decompNext.M);
    
    params.initPts = boQs;
    params.initVals = boVals;
  end

end

end
