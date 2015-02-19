function [func, funcProperties] = get2ModalFunction(numDims)

  bounds = repmat([-1, 1], numDims, 1);
  funcProperties.bounds = bounds;
  funcProperties.gaussVar = 0.01 * numDims^0.1;
  funcProperties.covarDiag = funcProperties.gaussVar * ones(1, numDims);
  funcProperties.centres12 = [0.72 -0.28; 0.38 0.18];
  funcProperties.centresRest = [-0.66 -0.18]';
  funcProperties.centres = ...
    [funcProperties.centres12, repmat(funcProperties.centresRest, 1, numDims -2)];
  probs = [0.3; 0.7];
%   probs = [0.3; 0.5; 0.3];
  funcProperties.centreProbs = probs;
  func = @(t) max( -700, ...
    log( ...
    probs(1) * mvnpdf(t, funcProperties.centres(1, :), funcProperties.covarDiag ) +...
    probs(2) * mvnpdf(t, funcProperties.centres(2, :), funcProperties.covarDiag ) ) ...
    );
  [~, maxCentreIdx] = max(probs);
  funcProperties.maxPt = funcProperties.centres(maxCentreIdx, :);
  funcProperties.maxVal = func(funcProperties.maxPt);

end

