function [func, funcProperties] = getAdditiveFunction(numDims, numDimsPerGroup)
% Returns an Additive Function that could be optimized. 

  if ~exist('numDimsPerGroup', 'var')
    numDimsPerGroup = min(5, round(numDims/2));
  end

  d = numDimsPerGroup;
  numRemDims = mod(numDims, d);
  D = numDims - numRemDims;
  numGroups = D/d;

  % Use this as a sub function for the three coordinates
  if d == 1
    [subFunc, subProperties] = get2Modal1DFunction;
  else
    [subFunc, subProperties] = get3ModalFunction(d);
  end

  func = @(X) objFunction(X, subFunc, D, d);
  % The bounds of the function for the D/d groups will be taken from subFunc.
  % The rest are [-1, 1]
  funcProperties.bounds = [ repmat( subProperties.bounds, numGroups, 1); ...
                            repmat( [-1 1], numRemDims, 1) ];
  funcProperties.maxPt = [ repmat( subProperties.maxPt, 1, numGroups), ...
                           repmat( [0], 1, numRemDims) ];
  funcProperties.maxVal = func(funcProperties.maxPt);
  funcProperties.decomposition = {};
  for k = 1:numGroups
    funcProperties.decomposition{k} = ((k-1)*d+1) : (k*d);
  end

end


% This is the actual function that does the computation
function val = objFunction(X, subFunc, D, d)
  X = X(:, 1:D); % discard the remaining dimensions
  numGroups = D/d;
  val = 0;
  for k = 1:numGroups
    val = val + subFunc( X(:, (k-1)*d+1: k*d) ); 
  end
end

