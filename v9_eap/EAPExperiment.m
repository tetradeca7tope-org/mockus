classdef EAPExperiment < handle

  properties
    allAntennas = { ...
      'A1 - 75mono fixed.nec', ...
      'A2 - 150mono fixed.nec', ...
      'A3 - 150dipole fixed.nec', ...
      'A4-300mono-fixed.nec', ...
      'A5 - 300dipole fixed.nec'};
    allSettings = { ...
      'S1 - 4ft square plate fixed.nec', ...
      'S2 - 4ft square plate w hole fixed.nec', ...
      'S3 - 4ft square plate w box fixed.nec', ...
      'S4 - 4ft square plate w sides fixed.nec', ...
      'S5 - 4ft square plate w box and sides fixed.nec', ...
      'S6 - 4ft square plate w box and sides sloped front fixed.nec', ...
      'S7 - 4ft square plate w sides sloped front fixed.nec'};

    alpha = 1;
    beta = 1;
    numAntennas;
    antennas;
    numDims;
    antennaSpaceBounds;
    problemSpaceBounds;
    settingIdx;
    setting;
    bounds;
    % EAP library details
%     eapDir = '../../eap/';
    eapDir = './sim/';
    luasDir;
    resultsDir;
    runsDir;
    mcFile;
    binName = 'evol_ant';

    % For executing the simulator
    paramsStr = ...
      ['params = {mutation = 0.0,exp_weight = 2,algorithm = "EX",', ...
       'run_simulator = 1,max_gain = 1,max_coup = 1,min_coup = 1,', ...
       'auto_seed = 1}'];
  end

  methods

    % Constructor
    function obj = EAPExperiment(settingIdx, numAntennas)

      % prelims
      obj.antennaSpaceBounds = repmat([-1, 1], 3*numel(obj.allAntennas), 1);
      % EAP stuff
%       obj.luasDir = sprintf('%sluas/', obj.eapDir);
      obj.luasDir = 'luas/';
      obj.resultsDir = sprintf('%sresults/', obj.eapDir);
      obj.runsDir = sprintf('%sruns/', obj.eapDir);

      % Set based on number of antennas
      if ~exist('settingIdx', 'var')
        settingIdx = 5;
      end
      obj.settingIdx = settingIdx;
      obj.setting = obj.allSettings{settingIdx};
      if ~exist('numAntennas', 'var')
        numAntennas = numel(obj.allAntennas);
      end
      obj.numAntennas = numAntennas;
      obj.antennas = obj.allAntennas(1:numAntennas);
      obj.numDims = 3 * numAntennas;
      obj.problemSpaceBounds = obj.antennaSpaceBounds(1:3*numAntennas, :);
      obj.bounds = repmat([0 1], obj.numDims);
      obj.mcFile = sprintf('%sind000000000a%02d.out', obj.runsDir, numAntennas);

    end

    function normCoords = getNormCoords(obj, trueCoords)
      normCoords = getNormParams(trueCoords, obj.problemSpaceBounds);
    end

    function trueCoords = getTrueCoords(obj, normCoords)
      trueCoords = getUnNormParams(normCoords, obj.problemSpaceBounds);
    end

    function [negFitnessVals] = normCoordFitness(obj, evalPts)
      negFitnessVals = obj.trueCoordFitness(obj.getTrueCoords(evalPts));
    end

    function [negFitnessVals] = trueCoordFitness(obj, evalPts)

      numEvalPts = size(evalPts, 1);
      negFitnessVals = zeros(numEvalPts, 1);

      for evalPtIter = 1:numEvalPts

        % Execute simulator here
        evalPt = evalPts(evalPtIter, :);
        success = obj.runSimulator(evalPt);
        if ~success
          error('Could not execute binary.\n');
        end

        % Read results
        currMCVal = obj.mutualCoupling();
        currGPVal = obj.gainPattern();
        negFitnessVals(evalPtIter) = -obj.alpha*currMCVal - obj.beta*currGPVal;
      end

    end % end trueCoordFitness


    % Run Simulator
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function success = runSimulator(obj, evalPt)

      % First write antenna locations to file
      luaFileStr = sprintf('load_platform(\"input/%s\")\n', obj.setting);
      for i = 1:obj.numAntennas
        currStr = sprintf(['add_antenna(\"input/%s\")\n', ...
          'add_point( %0.7f, %0.7f, %0.7f )\n'], ...
          obj.antennas{i}, ...
          evalPt(3*(i-1) + 1), evalPt(3*(i-1) + 2), evalPt(3*(i-1) + 3) );
        luaFileStr = sprintf('%s%s', luaFileStr, currStr);
      end
      luaFileStr = sprintf('%s%s', luaFileStr, obj.paramsStr);
      inFileStr = sprintf('%s%sin.lua', obj.eapDir, obj.luasDir);
      inFile = fopen(inFileStr, 'w');
      fprintf(inFile, luaFileStr);

      % Now execute the simulator
      commandStr = sprintf(...
        ['export LD_LIBRARY_PATH="/usr/lib/gcc/x86_64-linux-gnu/4.8" ; ', ...
        'cd sim && ./%s -i %sin.lua && cd ..'], ...
        obj.binName, obj.luasDir);
      commandStr,
      success = ~system(commandStr);

      success = 1;
    end

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

      % Pass the next two lines too
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

