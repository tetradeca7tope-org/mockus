% plot VJ results

collectVJResults;

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
experimentIdxs = boAddMaxVals(:, end, 1) > 80; % Discard the bad stuff
numExperiments = sum(experimentIdxs);
numExperiments,

numDimsPerGroupCands,
candPlotIdxs = [1 2 5 7 8 9]; %[2, 3, 4, 5]; % for 10D example
boAddMaxVals = boAddMaxVals(:,:,candPlotIdxs);
numDimsPerGroupCands = numDimsPerGroupCands(candPlotIdxs);

% First remove the zero entries
boAddMaxVals = boAddMaxVals(experimentIdxs, :, :);
randMaxVals = randMaxVals(experimentIdxs, :);

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
    figTitlePrefix = 'Max Vals';

  else % Now do Cumulative Regret

  end

  % Correct extra terms for diRect
  diRectReg = diRectReg(1:totalNumQueries);
  
  % Store the minimum and maximum for Plotting
  addRegMeanMinVals = AddRegMean(1, 1, :);
  minPlotVal = min([randRegMean(1); addRegMeanMinVals(:)]);
  addRegMeanMaxVals = AddRegMean(1, end, :);
  maxPlotVal = max([ randRegMean(end);addRegMeanMaxVals(:)]);

  % First plot the nominal value
  figure;
  plotFunc([0 totalNumQueries], [88.4 88.4], 'b--');
  hold on,

  % Plot Iteration statistics
  plotFunc(qqq, randRegMean(qqq), plotShapes{3}, 'Color', plotColours{3}, ...
     'MarkerSize', MARKER_SIZE, 'LineWidth', LINE_WIDTH); hold on,
  if regIter == 1, % Don't plot Cum Regret for DiRect
    plotFunc(qqq, diRectReg(qqq), plotShapes{4}, 'Color', plotColours{4}, ...
      'MarkerSize', MARKER_SIZE, 'LineWidth', LINE_WIDTH); hold on,
  end
  if regIter == 1
    legEntries = {'OpenCV', 'Random', 'DiRect' };
  else
    legEntries = {'Random', 'BO-KD', 'BO-UD'};
  end
  numBaseLegEntries = numel(legEntries);
  for i = 1:numel(candPlotIdxs)
    plotFunc(qqq, AddRegMean(1,qqq,i), plotShapes{4+i}, 'Color', plotColours{4+i}, ...
  'MarkerSize', MARKER_SIZE, 'LineWidth', LINE_WIDTH);
    numGroups = ceil(numDims/numDimsPerGroupCands(i));
    if i==1
    legEntries{numBaseLegEntries+i} = 'GP-UCB';
    else
    legEntries{numBaseLegEntries+i} = sprintf('Add-%d/%d', ...
        numDimsPerGroupCands(i), numGroups);
    end
  end
  legend(legEntries);
  % Now reproduce the curve without the bullets
  plotFunc(qq, randRegMean, 'Color', plotColours{3}, 'LineWidth', LINE_WIDTH); hold on,
  if regIter == 1, % Don't plot Cum Regret for DiRect
    plotFunc(qq, diRectReg, 'Color', plotColours{4}, 'LineWidth', LINE_WIDTH); hold on,
  end
  for i = 1:numel(candPlotIdxs)
    plotFunc(qq, AddRegMean(1,:,i), 'Color', plotColours{4+i}, 'LineWidth', LINE_WIDTH);
  end

  % Plot Error Bars
  if PLOT_ERR_BARS & (numExperiments > 1)
    errorbar(qqq, randRegMean(qqq), randRegStdErr(qqq), '.', 'Color', plotColours{3});
    for i = 1:numel(candPlotIdxs)
      errorbar(qqq, AddRegMean(1,qqq,i), AddRegStdErr(1,qqq,i), '.', ...
        'Color', plotColours{4+i});
    end
  end

  plotRange = maxPlotVal - minPlotVal; 
  xlim([0 1.05*totalNumQueries]);
  ylim([minPlotVal - plotRange*0.05, maxPlotVal + plotRange*0.05]);
  ylabel('Classification Accuracy');
  xlabel('Number of Queries (T)');
%   titlestr = 'Viola & Jones'; title(titlestr);

  set(0,'defaultAxesFontName', 'Dejavu Sans')
    set(findall(gca, '-property', 'FontSize'), 'FontSize', 18, ...
      'fontWeight', 'bold');

end

