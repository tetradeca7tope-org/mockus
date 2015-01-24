classdef LRGExperiment < handle

  properties
% %   Worked !
%     problemSpaceBounds = [-.01    .03 ; ... Omega_k: 0
%                           0.75     0.8   ; ... Omega_Lambda: 0.762
%                           0.10  0.105; ... omega_c: 0.1045
%                           0.02 0.03 ; ... Omega_b: 0.02233
%                           0.95   0.96 ; ... n_s: 0.951
%                           0.65  0.7; ... A_s: 0.6845
%                           -0.02  0.01 ; ... alpha: 0
%                           1.9     1.91 ; ...  b: 1.908
%                           30 31]; % QNL: 30.81
    problemSpaceBounds = [-.01    .03 ; ... Omega_k: 0
                          0.75     0.8   ; ... Omega_Lambda: 0.762
                          0.10  0.105; ... omega_c: 0.1045
                          0.02 0.03 ; ... Omega_b: 0.02233
                          0.95   0.96 ; ... n_s: 0.951
                          0.65  0.7; ... A_s: 0.6845
                          -0.02  0.01 ; ... alpha: 0
                          1.0     2.0 ; ...  b: 1.908
                          30 31]; % QNL: 30.81
%   Original
%     problemSpaceBounds = [-1    0.9 ; ... Omega_k: 0
%                           0     1   ; ... Omega_Lambda: 0.762
%                           0     1.2 ; ... omega_c: 0.1045
%                           0.001 0.25 ; ... Omega_b: 0.02233
%                           0.5   1.7 ; ... n_s: 0.951
%                           0.65  0.75; ... A_s: 0.6845
%                           -0.1  0.1 ; ... alpha: 0
%                           1.9   1.91  ; ...  b: 1.908
%                           30 31]; % QNL: 30.81
    % Just so that it doesn't favour diRect
    irrelOptPt = 0.1;
    irrelPtPenalty = 60;
    numDims;
    bounds;
    coordOrder;
    invertOrder;
  end

  methods

    % Constructor
    function obj = LRGExperiment(numDims, coordOrder)
      obj.numDims = numDims;
      obj.bounds = [repmat([0 1], numDims-9, 1); obj.problemSpaceBounds];
      obj.bounds = obj.bounds(coordOrder, :);
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

