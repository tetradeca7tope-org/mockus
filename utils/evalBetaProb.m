function [prob_vals] = evalBetaProb(a, b, x)
% evaluates the value of the pdf of a Beta(a, b) distribution at points x.

  prob_vals = x.^(a-1) .* (1-x).^(b-1) / beta(a, b);
end
