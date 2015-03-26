function [score, acc] = vjWrap(XTrain, YTrain, HaarCascade, params, bounds,...
options, sampleSize)

  if ~isfield(options, 'numTrain')
    numTrain = size(YTrain, 1);
  else
    numTrain = options.numTrain;
  end
  if ~isfield(options, 'doubleParams')
    doubleParams = false;
  else
    doubleParams = options.doubleParams;
  end

  XTrain = XTrain(:,:, 1:numTrain);
  YTrain = YTrain(1:numTrain);
  n = size(params,1);

  if doubleParams
    numParams = 44;
  else
    numParams = 22;
  end
  
  assert(numParams == size(params,2), 'dimension do not match'); 
  thresholds = getUnNormParams(params, bounds);

  %%%%%%%%% use the best lower bounds 
  %bestParams = [0.8227 6.9566 9.4985 18.4130 15.3241 21.0106 23.9188 24.5279 ... 
  %27.1534 34.5541 39.1073 50.6105 54.6201 50.1697 66.6691 67.6989 69.2288 ...
  %79.2491 87.6960 90.2533 104.7492 105.7611];   
  %best = [bestParams bestParams.*2];
  %thresholds = best;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  score = zeros(n, 1);
  acc = zeros(n, 1);
  for i = 1:n
    [s, a] = trainThresholds(XTrain, YTrain, HaarCascade, thresholds(i,:)');
    score(i) = s;
    acc(i) = a;
  end
end
