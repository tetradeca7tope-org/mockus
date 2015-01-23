% v16: SDSS

clear all;
close all;
rng('shuffle');
addpath ../GPLibkky/
addpath ../UCLibkky/
addpath ../LipLibkky/
addpath ../MCMCLibkky/
addpath ../ABCLibkky/
addpath ../helper/

warning off;

% Load all constants
loadConstants;

% Don't store any results
% numResultsToBeStored = 0;
% numMCMCResultsToBeStored = 0;

for experimentIter = 1:NUM_EXPERIMENTS

  fprintf('\nExperiment #%d\n', experimentIter);

  fprintf('Uncertainty Reduction\n');
  UncertaintyReductionForLRG;

%   fprintf('Max Band Point\n');
%   MaxBandPointForLRG;

  fprintf('MCMC\n');
  MCMCForLRG;

  fprintf('RAND\n');
  RANDForLRG;

  fprintf('\n');

  save( saveFileName, 'uc_errs', 'mbp_errs', 'rand_errs', 'mcmc_errs', ...
    'mcmcReg_errs');
  savePtsFileName = sprintf('%s_%d.mat', savePtsFilePrefix, experimentIter);
  save( savePtsFileName, 'ucPts', 'ucLogProbs', 'mcmcQueries', ...
  'mcmcLogProbs', 'mrPts', 'mrLogProbs', 'randPts', 'randLogProbs', ...
  'mcmcSamples');
end


gtSaveFileName = sprintf('gt10_%s.mat', timestamp);
save(gtSaveFileName, 'ucPts', 'ucLogProbs', 'mcmcQueries', ...
  'mcmcLogProbs', 'mrPts', 'mrLogProbs', 'randPts', 'randLogProbs', ...
  'mcmcSamples');
system('beep');
