% Generates some figures for BO/GPB
close all;
clear all;
addpath ~/libs/kky-matlab/GPLibkky/
addpath ~/libs/kky-matlab/utils/
rng(104);

N = 100;
numSamples = 20;

PLOT_UCB_ON_FUNC = true;
% PLOT_UCB_ON_FUNC = false;

func = @(t) - 70*(t-0).* (t-0.35).* (t+0.55).* (t-0.65).* (t-0.98) - 0.27;

numPts = 25;
stopPts = [1:15, numPts];

numPts = 104;
% stopPts = [1:10:numPts];

numCandPts = 300;
candPts = linspace(0,1,numCandPts+2); candPts = candPts(2: (end-1))';
gpConfWidth = 3;

LW = 3;
MS = 25;
LWS = 1.5;
FS = 26;

plotlims = [0 1 -0.75 2.4];
% figSize = [85 261 898 420];
figSize = [85 061 898 500];
saveFormat = 'eps';
yLabPosn = [-0.05 1.90];
xLabPosn = [0.95 -0.80];

c1 = [150 75 0]/255;
c2 = 'r';

% saveas(gcf, 'func', saveFormat);
% ylabh = get(gca, 'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') - [0.2 0])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine Kernel parameters
% bw = 0.11; scale = 0.18509; mu0 = 0; noise = 0.001; %0.05;
bw = 0.12; scale = 0.18509; mu0 = 0; noise = 0.001; %0.05;

% Initialisation
muCandPts = mu0 * ones(numCandPts, 1);
Kte = scale * GaussKernel(bw, candPts, candPts);
sigmaCandPts = sqrt(diag(Kte));
betath = 2;
currEvalPts = zeros(0,1);
currEvalVals = zeros(0,1);

% Plot the original GP
    top = muCandPts + gpConfWidth*sigmaCandPts;
    bottom = muCandPts - gpConfWidth*sigmaCandPts;
    A = [candPts' fliplr(candPts')]; B = [top' fliplr(bottom')]; h = fill(A, B, [0.9 0.9 0.9]);
    set(h, 'EdgeColor', 'None'); hold on,
    %  NOw plot the function
    th = linspace(0,1,100)';
    [maxVal, maxIdx] = max(func(th));
    maxPt = th(maxIdx);
    plot(th, func(th), 'LineWidth', LW, 'Color', 'k'); hold on,
    set(gca, 'Ytick', []);
    set(gca, 'Xtick', []);
    axis(plotlims);
    xlabel('$x$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
    yLabPosn = [-0.05 0.8*plotlims(4)];
    ylabel('$f(x)$', 'rot', 0, 'Position', yLabPosn, 'FontSize', FS, ...
      'Interpreter', 'Latex');
    set(gcf, 'Position', figSize);
%     pause,
    set(gca, 'Position', [0.08 0.08 0.91 0.91], 'units', 'normalized')
    box off,
    box on,
    pause

for t = 1:numPts

  % Will need this to plot the uncertainty region
  prevMuCandPts = muCandPts;
  prevSigmaCandPts = sigmaCandPts;

  betat = 5*log(t + exp(1));
%   betath = 2*log(1 + exp(1));
  betath = sqrt(betat);
  ucbCandPts = muCandPts + betath * sigmaCandPts;
  if t == 1,
    nextPtUCB = mu0 + betath * sigmaCandPts(1);
    nextPt = 0.5;
  else
    [nextPtUCB, nextPtIdx] = max(ucbCandPts);
    nextPt = candPts(nextPtIdx);
  end


  % Evaluate function GP
  nextVal = func(nextPt);
  currEvalPts = [currEvalPts; nextPt];
  currEvalVals = [currEvalVals; nextVal];

  % Update GP
  Ktr = scale * GaussKernel(bw, currEvalPts) + noise * eye(t);
  L = stableCholesky(Ktr);
  alpha = L' \ (L \ (currEvalVals - mu0 * ones(t,1) ));
  Ktetr = scale * GaussKernel(bw, candPts, currEvalPts);
  muCandPts = mu0 + Ktetr * alpha;
  V = L \ Ktetr';
  KCandPts = Kte - V'*V;
  sigmaCandPts = sqrt(diag(KCandPts));

  if any(stopPts == t)

    if ~PLOT_UCB_ON_FUNC
      % Plot the figures
      figure;
      if t == 1
        modUCBVals = 0.5 * ones(numCandPts, 1);
        modNextPtUCB = 0.5;
      else
        modUCBVals = 0.05 + ...
          0.88 * (ucbCandPts - min(ucbCandPts))/(max(ucbCandPts) - min(ucbCandPts));
        modNextPtUCB = modUCBVals(nextPtIdx);
      end
      plot(candPts, modUCBVals, 'Color', 'b', 'LineWidth', LW);  hold on,
      plot(nextPt, modNextPtUCB, 'r*', 'LineWidth', 1.5*LW, 'MarkerSize', MS);

      set(gca, 'Xtick', []);
      set(gca, 'Ytick', []);
      axis([0 1 0 1]);
      % set(gca, 'position', [0 0 1 1],'units','normalized');
      xlabel('$x$', 'Position', [0.95, -0.04], 'FontSize', FS, 'Interpreter','Latex');
      ylabel('$\varphi_t(x)$', 'Interpreter', 'Latex', ...
         'rot', 0, 'Position', [-0.06 0.9], 'FontSize', FS);
      set(gcf, 'Position', figSize);
    end


    % The function itself
    figure;
%     top = muCandPts + gpConfWidth*sigmaCandPts;
%     bottom = muCandPts - gpConfWidth*sigmaCandPts;
    top = prevMuCandPts + gpConfWidth*prevSigmaCandPts;
    bottom = prevMuCandPts - gpConfWidth*prevSigmaCandPts;
    A = [candPts' fliplr(candPts')];
    B = [top' fliplr(bottom')];
    h = fill(A, B, [0.9 0.9 0.9]);
%     h = fill(A, B, [173 255 47]/255);
    set(h, 'EdgeColor', 'None'); hold on,

    plot(th, func(th), 'LineWidth', LW, 'Color', 'k'); hold on,
    plot(currEvalPts(1:end-1), currEvalVals(1:end-1), 'x', 'Color', 'm', ...
      'LineWidth', LW, 'MarkerSize', 0.7*MS);

    % Plot the UCB
    if PLOT_UCB_ON_FUNC
      plot(candPts, ucbCandPts, '--', 'Color', 'b', 'LineWidth', LW);  hold on,
      plot(nextPt, nextPtUCB, '*', 'Color', [255, 100, 0]/255, ...
        'LineWidth', 1.5*LW, 'MarkerSize', 0.8*MS);
      plot([nextPt, nextPt], [-0.75, nextPtUCB], '--', ...
        'Color', [255, 100, 0]/255, ...
        'LineWidth', 0.8*LW, 'MarkerSize', 0.8*MS);
    end

%     if t ~= numPts
%       plot(nextPt, nextVal, 'r*', 'LineWidth', 2*LW, 'MarkerSize', MS);
%     end

    timeStr = sprintf('$t=%d$', t);
    text(0.05, -0.6, timeStr, 'FontSize', FS, ...
      'Interpreter', 'Latex');
    set(gca, 'Ytick', []);
    set(gca, 'Xtick', []);
%     if t <= 30,
%       plotlims = [0 1 -0.75 2.4];
% %       yLabPosn = [-0.05 1.90];
%     else
%       plotlims = [0 1 -0.75 1.2];
% %       yLabPosn = [-0.05 0.8];
%     end

    yLabPosn = [-0.05 0.8*plotlims(4)];

    axis(plotlims);
    xlabel('$x$', 'Position', xLabPosn, 'FontSize', FS, 'Interpreter', 'Latex');
    ylabel('$f(x)$', 'rot', 0, 'Position', yLabPosn, 'FontSize', FS, ...
      'Interpreter', 'Latex');
    set(gcf, 'Position', figSize);
    set(gca, 'Position', [0.08 0.08 0.91 0.91], 'units', 'normalized')
    box off;
    box on;



%     pause; close all;
  end
  


end


