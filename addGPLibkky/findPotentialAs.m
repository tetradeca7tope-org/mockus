function [As] = findPotentialAs(func, D, p)
% Returns a set of numAs candidates for the A matrix.
% We will randomly initialize at numInitPts points and then regroup the
% variables.

  numAs = 10; 
  numInitPts = min(1000, 10*d*p);


  % Set Steifel Optimization parameters
  stiefelOpts.record = 0;
  stiefelOpts.xtol = 1e-5;
  stiefelOpts.gtol = 1e-5;
  sti

  % TODO: See if the following strategy is useful: locally optimize, then
  % project to the set of permutation matrices, then locally optimize again,
  % project and repeat so on...

end
