function A = decompOptPartial(func, D, d, M)

  numTrials = D*d*M;
  currBestVal = inf;

  A = getRandPermMat(D);
  for i = 1:numTrials
    P = getRandPermMat(D);
    val = func(P);
    if val < currBestVal
      A = P;
      currBestVal = val;
    end
  end

end

