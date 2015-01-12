function [simpleRegret, cumRegretByT, cumRegret] = getRegrets(maxVal, histories)

  numQueries = size(histories, 1);

  diffs = maxVal - histories;
  cumRegret = cumsum(diffs);
  cumRegretByT = cumRegret ./ (1:numQueries)';

  simpleRegret = zeros(numQueries, 1);
  simpleRegret(1) = diffs(1);
  for i = 2:numQueries
    simpleRegret(i) = min( simpleRegret(i-1), diffs(i) );
  end

end

