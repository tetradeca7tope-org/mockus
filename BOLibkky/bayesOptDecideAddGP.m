function [maxVal, maxPt, boQueries, boVals, history] = ...
bayesOptDecideAddGP(oracle, decomp, bounds, numIters, params)

% decomp is a cell array, each cell is a struct with the following 
% fields: d, M, group, score. 

DECOMP_KNOWN = 'known';
DECOMP_LEARN = 'learn';
DECOMP_RAND = 'random';
DECOMP_PLEARN = 'partialLearn';
DECOMP_DECIDE = 'decide';

if ~strcmp(params.decompStrategy, DECOMP_DECIDE)
  [maxVal, maxPt, boQueries, boVals, history] = ...
  bayesOptDecompAddGP(oracle, decomp, bounds, numIters, params);

else
  boQueries = [];
  boVals = [];
  history = [];

  numDecompChoices = numel(decomp);
  numInnerIters = floor(numIters / numDecompChoices);

  for i=1:numDecompChoices
    thisDecomp.d = decomp{i}.d;
    thisDecomp.M = decomp{i}.M;
    params.decompStrategy = DECOMP_PLEARN;

    [maxVal, maxPt, thisBoQueries, thisBoVals, thisHistory]=...
    bayesOptDecompAddGP(oracle, thisDecomp, bounds, numInnerIters, params);
   
    size(thisBoQueries);
    size(thisBoVals);
    if (i==1)
      boQueries = [boQueries;thisBoQueries];
      boVals = [boVals;thisBoVals];
    else 
      histLen = numel(thisBoVals);
      boQueries = [boQueries;thisBoQueries(histLen-numInnerIters+1:end,:)];
      boVals = [boVals;thisBoVals(histLen-numInnerIters+1:end)];
    end

    params.initPts = boQueries;
    params.initVals = boVals; 
  end
  
  size(boVals)
end
end