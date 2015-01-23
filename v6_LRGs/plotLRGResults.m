% plot VJ results

% saveFileName = 
load(saveFileName);
% Prelims
close all;
PLOT_ERR_BARS = true;
% PLOT_ERR_BARS = false;
NUM_ERR_BARS = 10;
MARKER_SIZE = 8;
LINE_WIDTH = 2;
SAVE_FILE_FORMAT = 'png';

resultsDir = 'results/';
plotColours = {'c', 'b', 'r', 'm', 'k', 'g', [255 128 0]/255, ...
  [76, 0, 153]/253, [102 102 0]/255, 'c', [76, 160, 163]/253, 'b', 'y'};
plotShapesDot = {'o.', '+.', '*.', 'x.', 's.', 'd.', '^.', 'p.', '>.', 'v.', ...
      'o.', '+.', '*.', 'x.', 's.'};
plotShapes = {'o', '+', '*', 'x', 's', 'd', '^', 'p', '>', 'v', ...
      'o', '+', '*', 'x', 's'};
plotFunc = @semilogx;
plotFunc = @loglog;
plotFunc = @semilogy;
% plotFunc = @plot
qq = 1:totalNumQueries;
qqq = round(linspace(1,totalNumQueries, NUM_ERR_BARS+2)); qqq = qqq(2:end-1);
numExperiments = sum((randMaxVals(:,1) ~= 0));
numExperiments,


% First remove the zero entries
boAddMaxVals = boAddMaxVals(1:numExperiments, :, :);
randMaxVals = randMaxVals(1:numExperiments, :);
% Cum Regrets
boAddCumRewards = boAddCumRewards(1:numExperiments, :, :);
randCumRewards = randCumRewards(1:numExperiments, :);

for regIter = 1:1

  if regIter == 1 % First Do Simple Regret
    % Mean
    AddRegMean = mean(boAddMaxVals, 1);
    randRegMean = mean(randMaxVals, 1);
    % Std
    AddRegStdErr = std(boAddMaxVals, 1)/sqrt(numExperiments);
    randRegStdErr = std(randMaxVals, 1)/sqrt(numExperiments);
    % For diRect
    diRectReg = diRectMaxVals;
%     diRectTail = 400; diRectReg(diRectTail:end) = diRectReg(diRectTail);
    figTitlePrefix = 'Max Vals';

%     buggyPts = UDRegMean < KDRegMean;
%     KDRegMean(buggyPts) = 0.5*KDRegMean(buggyPts) + 0.5*UDRegMean(buggyPts);

  else % Now do Cumulative Regret
    % Mean
    KDRegMean = mean(boKDCumRegrets, 1);
    UDRegMean = mean(boUDCumRegrets, 1);
    AddRegMean = mean(boAddCumRegrets, 1);
    randRegMean = mean(randCumRegrets, 1);
    % Std
    KDRegStdErr = std(boKDCumRegrets, 1)/sqrt(numExperiments);
    UDRegStdErr = std(boUDCumRegrets, 1)/sqrt(numExperiments);
    AddRegStdErr = std(boAddCumRegrets, 1)/sqrt(numExperiments);
    randRegStdErr = std(randCumRegrets, 1)/sqrt(numExperiments);
    % For diRect
%     diRectReg = diRectCumRegret;
    figTitlePrefix = 'Cumulative-Regret';

%     buggyPts = UDRegMean < KDRegMean;
%     KDRegMean(buggyPts) = 0.5*KDRegMean(buggyPts) + 0.5*UDRegMean(buggyPts);

  end

  % strip diRect
  diRectReg = diRectReg(1:totalNumQueries, :);

  % Plot Iteration statistics
  figure;
  plotFunc(qqq, randRegMean(qqq), plotShapes{3}, 'Color', plotColours{3}, ...
     'MarkerSize', MARKER_SIZE, 'LineWidth', LINE_WIDTH); hold on,
  if regIter == 1, % Don't plot Cum Regret for DiRect
    plotFunc(qqq, diRectReg(qqq), plotShapes{4}, 'Color', plotColours{4}, ...
      'MarkerSize', MARKER_SIZE, 'LineWidth', LINE_WIDTH); hold on,
  end
  if regIter == 1
    legEntries = {'Random', 'DiRect' };
  else
    legEntries = {'Random', 'BO-KD', 'BO-UD'};
  end
  numBaseLegEntries = numel(legEntries);
  for i = 1:numdCands
    plotFunc(qqq, AddRegMean(1,qqq,i), plotShapes{4+i}, 'Color', plotColours{4+i}, ...
  'MarkerSize', MARKER_SIZE, 'LineWidth', LINE_WIDTH);
    legEntries{numBaseLegEntries+i} = sprintf('BO-%d', numDimsPerGroupCands(i));
  end
  legend(legEntries);
  % Now reproduce the curve without the bullets
  plotFunc(qq, randRegMean, 'Color', plotColours{3}, 'LineWidth', LINE_WIDTH); hold on,
  if regIter == 1, % Don't plot Cum Regret for DiRect
    plotFunc(qq, diRectReg, 'Color', plotColours{4}, 'LineWidth', LINE_WIDTH); hold on,
  end
  for i = 1:numdCands
    plotFunc(qq, AddRegMean(1,:,i), 'Color', plotColours{4+i}, 'LineWidth', LINE_WIDTH);
  end

  % Plot Error Bars
  if PLOT_ERR_BARS & (numExperiments > 1)
    errorbar(qqq, randRegMean(qqq), randRegStdErr(qqq), '.', 'Color', plotColours{3});
    for i = 1:numdCands
      errorbar(qqq, AddRegMean(1,qqq,i), AddRegStdErr(1,qqq,i), '.', ...
        'Color', plotColours{4+i});
    end
  end

  xlim([0 1.1*totalNumQueries]);
%   ylim([60 100]);
  titleStr = sprintf('%s, D = %d', figTitlePrefix, numDims);
  title(titleStr);

  % Also Save the image
  imFileName = sprintf('%s-%s.%s', saveFileName(1:(end-4)), ...
    figTitlePrefix, SAVE_FILE_FORMAT);
  saveas(gcf, imFileName, SAVE_FILE_FORMAT);

end
