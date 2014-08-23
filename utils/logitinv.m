function y = logitinv(x)
  y = exp(x) ./ (1 + exp(x));
end
