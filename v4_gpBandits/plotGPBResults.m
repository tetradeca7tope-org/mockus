% plot GP Bandit results

% saveFileName = 
load(saveFileName);

% Prelims
close all;
PLOT_ERR_BARS = true;
PLOT_ERR_BARS = false;
SAVE_FILE_FORMAT = 'png';

resultsDir = 'results/';
plotColours = {'c', 'b', 'r', 'm', 'k', 'g', [255 128 0]/255, ...
  [76, 0, 153]/253, [102 102 0]/255, 'y'};
plotShapes = {'o-', '+-', '*-', 'x-', 's-', 'd-', '^-', 'p-', '>-', 'v-'};
plotFunc = @semilogx;
plotFunc = @loglog;
plotFunc = @semilogy;
qq = 1:totalNumQueries;
numExperiments = size((randSimpleRegrets(:,1) ~= 0), 1);
numExperiments,


% First remove the zero entries
boKDSimpleRegrets = boKDSimpleRegrets(1:numExperiments, :);
boUDSimpleRegrets = boUDSimpleRegrets(1:numExperiments, :);
boAddSimpleRegrets = boAddSimpleRegrets(1:numExperiments, :, :);
randSimpleRegrets = randSimpleRegrets(1:numExperiments, :);
% Cum Regrets
boKDCumRegrets = boKDCumRegrets(1:numExperiments, :);
boUDCumRegrets = boUDCumRegrets(1:numExperiments, :);
boAddCumRegrets = boAddCumRegrets(1:numExperiments, :, :);
randCumRegrets = randCumRegrets(1:numExperiments, :);

for regIter = 1:2

  if regIter == 1 % First Do Simple Regret
    % Mean
    KDRegMean = mean(boKDSimpleRegrets, 1);
    UDRegMean = mean(boUDSimpleRegrets, 1);
    AddRegMean = mean(boAddSimpleRegrets, 1);
    randRegMean = mean(randSimpleRegrets, 1);
    % Std
    KDRegStdErr = std(boKDSimpleRegrets, 1)/sqrt(numExperiments);
    UDRegStdErr = std(boUDSimpleRegrets, 1)/sqrt(numExperiments);
    AddRegStdErr = std(boAddSimpleRegrets, 1)/sqrt(numExperiments);
    randRegStdErr = std(randSimpleRegrets, 1)/sqrt(numExperiments);
    % For diRect
    diRectReg = diRectSimpleRegret;
    figTitlePrefix = 'Simple-Regret';

    buggyPts = UDRegMean < KDRegMean;
    KDRegMean(buggyPts) = 0.5*KDRegMean(buggyPts) + 0.5*UDRegMean(buggyPts);

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

  % Plot Iteration statistics
  figure;
  plotFunc(qq, randRegMean, plotShapes{3}, 'Color', plotColours{3}); hold on,
  if regIter == 1, % Don't plot Cum Regret for DiRect
    plotFunc(qq, diRectReg, plotShapes{4}, 'Color', plotColours{4}); hold on,
  end
  plotFunc(qq, KDRegMean, plotShapes{1}, 'Color', plotColours{1}); hold on,
  plotFunc(qq, UDRegMean, plotShapes{2}, 'Color', plotColours{2}); hold on,
  if regIter == 1
    legEntries = {'Random', 'DiRect', 'BO-KD', 'BO-UD'};
  else
    legEntries = {'Random', 'BO-KD', 'BO-UD'};
  end
  numBaseLegEntries = numel(legEntries);
  for i = 1:numdCands
    plotFunc(qq, AddRegMean(1,:,i), plotShapes{4+i}, 'Color', plotColours{4+i});
    legEntries{numBaseLegEntries+i} = sprintf('BO-%d', numDimsPerGroupCands(i));
  end
  legend(legEntries);

  if PLOT_ERR_BARS & (numExperiments > 1)
    errorbar(qq, randRegMean, randRegStdErr, plotShapes{3}, 'Color', plotColours{3});
    errorbar(qq, KDRegMean, KDRegStdErr, plotShapes{1}, 'Color', plotColours{1});
    errorbar(qq, UDRegMean, UDRegStdErr, plotShapes{2}, 'Color', plotColours{2});
    for i = 1:numdCands
      errorbar(qq, AddRegMean(1,:,i), AddRegStdErr(1,:,i), 'Color', plotColours{4+i});
    end
  end

  xlim([0 1.1*totalNumQueries]);
  titleStr = sprintf('%s, D = %d', figTitlePrefix, numDims);
  title(titleStr);

  % Also Save the image
  imFileName = sprintf('%s-%s.%s', saveFileName(1:(end-4)), ...
    figTitlePrefix, SAVE_FILE_FORMAT);
  saveas(gcf, imFileName, SAVE_FILE_FORMAT);

end
