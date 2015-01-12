function [score, acc] = vjWrap(XTrain, YTrain, HaarCascade, params, options)

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
  n = size(params, 1);

% params is in [0,1]^22
  if doubleParams
    numParams = 44;
  else
    numParams = 22;
  end
  bounds = 110 * [zeros(numParams,1) ones(numParams,1)];

  thresholds = getUnNormParams(params, bounds);
  score = zeros(n, 1);
  acc = zeros(n, 1);
  for i = 1:n
    [s, a] = trainThresholds(XTrain, YTrain, HaarCascade, thresholds(i, :)');
    score(i) = s;
    acc(i) = a;
  end

end

