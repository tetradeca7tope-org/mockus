function [negNlml, nAG] = negNormMargLikelihood(sigmaSms, sigmaPrs, ...
  decomposition, A, X, y, meanFuncs, commonMeanFunc, noises, commonNoise) 

  if nargout == 2,
    [nlml, AG] = normRotMargLikelihood(sigmaSms, sigmaPrs, decomposition, A, ...
      X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
  else
    nlml = normRotMargLikelihood(sigmaSms, sigmaPrs, decomposition, A, ...
      X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
    AG = [];
  end
  negNlml = -nlml;
  nAG = -AG;

end

