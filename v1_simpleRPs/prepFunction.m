% prepares the function to be optimized

funcProperties.gaussVar = 0.5 * sqrt(numDims);
funcProperties.centres12 = [0.8 0.3; 0.6 0.8; 0.2 0.2];
funcProperties.centresRest = [0.2 0.5 0.8]';
funcProperties.centres = ...
  [funcProperties.centres12, repmat(funcProperties.centresRest, 1, numDims -2)];
func = @(t) log(
  0.2 * mvnpdf(t, funcProperties.centres(1, :), funcProperties.gaussVar) + ...
  0.5 * mvnpdf(t, funcProperties.centres(2, :), funcProperties.gaussVar) + ...
  0.3 * mvnpdf(t, funcProperties.centres(3, :), funcProperties.gaussVar) );

