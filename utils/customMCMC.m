function [samples, queries, logProbs] = ...
  CustomMCMC(numSamples, proposalStd, initPt, evalLogProb)
% Returns: the samples collected by MCMC, the pts at which the likelihood was
% quried and the logLikelihoods at those points.

  num_dims = size(initPt, 1);
  
  samples = zeros(numSamples, num_dims);
  queries = zeros(numSamples, num_dims);
  logProbs = zeros(numSamples, 1);

  currPt = initPt;
  currLogProb = evalLogProb(currPt');

  for sample_iter = 1:numSamples
%     if mod(sample_iter, 1000) == 0,
%       fprintf('MCMC Iter: %d\n', sample_iter);
%     end
    nextPt = currPt + proposalStd * randn(num_dims, 1);
    nextLogProb = evalLogProb(nextPt');
    queries(sample_iter, :) = nextPt;
    logProbs(sample_iter, :) = nextLogProb;

    log_prob_ratio = nextLogProb - currLogProb;
    if (log(rand()) < log_prob_ratio) % then accept the proposal
      currPt = nextPt;
      currLogProb = nextLogProb;
    end
    samples(sample_iter, :) = currPt';
  end

end
