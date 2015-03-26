function [score, acc] = trainThresholds(XTrain, YTrain, HaarCascade, ...
thresholds)

  numImages = size(YTrain, 1);   

  %allInfo = zeros(numImages, 44);
  %save('allInfo','allInfo');

  % Predicts
  predicts = zeros(numImages, 1);

  % Options
  Options.Resize = false;
  Options.Verbose = false;

  if isstr(HaarCascade)
    HaarCascade = GetHaarCasade(HaarCascade);
  end
    
  for i = 1:numImages
    objs = ObjectDetection( XTrain(:,:,i), HaarCascade, Options, thresholds,i);
    predicts(i) = ~isempty(objs); 
  end
  
  acc = sum(predicts == YTrain)/numImages;
  score = 100 * acc; 

end
