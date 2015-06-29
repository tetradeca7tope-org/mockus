classdef EAPExperiment < handle

  properties
    allAntennas = {'A1 - 750mono fixed.nec', ...
                'A2 - 150momo fixed.mec', ...
                'A3 - 150dipole fixed.nec', ...
                'A4-300mono-fixed.nec', ...
                'A5 - 300dipole fixed.nec'};

    alpha = 1;
    beta = 1;
    numAntennas;
    antennas;
    numDims;
    antennaSpaceBounds;
    problemSpaceBounds;
    bounds;
    % EAP library details
    eapDir = '../../eap/';
    luasDir;
    resultsDir;
    runsDir;
    mcFile;
  end

  methods

    % Constructor
    function obj = EAPExperiment(numAntennas)

      % prelims
      obj.antennaSpaceBounds = repmat([-1, 1], numel(obj.allAntennas), 1);
      % EAP stuff
      obj.luasDir = sprintf('%sluas/', obj.eapDir);
      obj.resultsDir = sprintf('%sresults/', obj.eapDir);
      obj.runsDir = sprintf('%sruns/', obj.eapDir);
      obj.mcFile = sprintf('%sind000000000a%02d.out', obj.runsDir, numAntennas);

      % Set based on number of antennas
      if ~exist('numAntennas', 'var')
        numAntennas = numel(obj.allAntennas);
      end
      obj.numAntennas = numAntennas;
      obj.antennas = obj.allAntennas(1:numAntennas);
      obj.numDims = 3 * numAntennas;
      obj.problemSpaceBounds = obj.antennaSpaceBounds(1:numAntennas, :);
      obj.bounds = repmat([0 1], obj.numDims);

    end

    function normCoords = getNormCoords(obj, trueCoords)
      normCoords = getNormParams(trueCoords, obj.bounds);
    end

    function trueCoords = getTrueCoords(obj, normCoords)
      trueCoords = getUnNormParams(normCoords, obj.bounds);
    end

    function [negFitnessVals] = normCoordFitness(obj, evalPts)
      negFitnessVals = obj.trueCoordLogJointProbs(obj.getTrueCoords(evalPts));
    end

    function [negFitnessVals] = trueCoordFitness(obj, evalPts)

      numEvalPts = size(evalPts, 1);
      negFitnessVals = zeros(numEvalPts, 1);

      for evalPtIter = 1:numEvalPts

        % Execute simulator here

        % Read results
        currMCVal = obj.mutualCoupling();
        currGPVal = obj.gainPattern();
        negFitnessVals(evalPtIter) = -obj.alpha*currMCVal - obj.beta*currGPVal;
      end

    end % end trueCoordFitness


    % Mutual Coupling
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function mcVal = mutualCoupling(obj)

      % Read file till you get to the MC values  
%       obj.mcFile = 'temp.txt';
      mcFileId = fopen(obj.mcFile, 'r');
      currIdWord = 'temp';
      while ~strcmp(currIdWord, 'COUPLING')
        currIdWord = obj.getIdWordInLine(mcFileId);
      end

      % Pass the next two files too
      fgets(mcFileId);
      fgets(mcFileId);

      mcVal = 0;
      for i = 1:nchoosek(obj.numAntennas, 2)
        q = textscan(fgets(mcFileId), '%s');
        mcVal = mcVal + str2num(q{1}{7});
      end

    end % end mutualCoupling

    % Gain Pattern
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function gpVal = gainPattern(obj)
      gpVal = 0; % return 0 for now.
    end

    % Ancillary functions
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Subroutine to read the first word of a line
    function idWord = getIdWordInLine(obj, fileId)
      currLine = fgets(fileId);
        words = textscan(currLine, '%s');
%       try 
%         words = textscan(currLine, '%s');
%       catch err
%         currLine,
%         pause;
%       end
      if numel(words{1})<2
        idWord = 'temp';
      else
        idWord = words{1}{2};
      end
    end



  end % methods
end % classdef

