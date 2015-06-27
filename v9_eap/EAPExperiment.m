classdef EAPExperiment < handle

  properties
    allAntennas = {'A1 - 750mono fixed.nec', ...
                'A2 - 150momo fixed.mec', ...
                'A3 - 150dipole fixed.nec', ...
                'A4-300mono-fixed.nec', ...
                'A5 - 300dipole fixed.nec'};
    antennaSpaceBounds = repmat([-1, 1], numel(allAntennas), 1);

    numAntennas;
    antennas;
    numDims;
    problemSpaceBounds;
    bounds;
  end

  methods

    % Constructor
    function obj = LRGExperiment(numAntennas)
      if ~exist('numAntennas', 'var')
        numAntennas = numel(antennas);
      end
      obj.numAntennas = numAntennas;
      obj.antennas = obj.allAntennas(1:numAntennas);
      obj.numDims = 3 * numAntennas;
      obj.problemSpaceBounds = obj.antennaSpaceBounds(1:numAntennas, :);
      obj.bounds = repmat([0 1], numDims);

%       obj.bounds = [repmat([0 1], numDims-9, 1); obj.problemSpaceBounds];
%       obj.bounds = obj.bounds(coordOrder, :);
    end

    function normCoords = getNormCoords(obj, trueCoords)
      normCoords = getNormParams(trueCoords, obj.bounds);
    end

    function trueCoords = getTrueCoords(obj, normCoords)
      trueCoords = getUnNormParams(normCoords, obj.bounds);
    end

    function [logJointProbs] = normCoordLogJointProbs(obj, evalPts)
      [logJointProbs] = obj.trueCoordLogJointProbs(obj.getTrueCoords(evalPts));
    end

    function [logJointProbs] = trueCoordLogJointProbs(obj, evalPts)
      % Identify points that are outside the domain
      evalPts = evalPts(:, obj.invertOrder);
      irrelPts = evalPts(:, 1:obj.numDims-9);
      modEvalPts = evalPts(:, obj.numDims-8:end);
      logJointProbs = lrgLogLiklWrap(modEvalPts);
      logJointProbs = logJointProbs - ...
        obj.irrelPtPenalty * sum((irrelPts - obj.irrelOptPt).^2, 2);
    end

  end % methods
end % classdef

