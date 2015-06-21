function [score, acc] = trainThresholds(XTrain, YTrain, HaarCascade, thresholds)

  numTrain = size(YTrain, 1);

  % Predicts
  predicts = zeros(numTrain, 1);

  % Options
  Options.Resize = false;
  Options.Verbose = false;

  if isstr(HaarCascade)
    HaarCascade = GetHaarCasade(HaarCascade);
  end
  
  for i = 1:numTrain
    objs = ObjectDetection( XTrain(:,:,i), HaarCascade, Options, thresholds);
    predicts(i) = ~isempty(objs); 
%     objs, predicts(i),
%     pause,
  end
  
  acc = sum(predicts == YTrain)/numTrain;

  % Modify it to obtain a smoother function ?
%   score = log(1 + 1e6*acc);
  score = 100 * acc; 

end

