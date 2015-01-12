% Unit test for vjWrap.m

load VJ/vjdata
addpath VJ/
addpath ../utils/
HaarCascade = GetHaarCasade('VJ/HaarCascades/haarcascade_frontalface_alt.mat');

bounds = 110*[zeros(22,1) ones(22,1)];
bestParams = [0.8227 6.9566 9.4985 18.4130 15.3241 21.0106 23.9188 24.5279 ... 
  27.1534 34.5541 39.1073 50.6105 54.6201 50.1697 66.6691 67.6989 69.2288 ...
  79.2491 87.6960 90.2533 104.7492 105.7611];
roundBestParams = round(bestParams);
% params = [  0.9*ones(1, 22); rand(1, 22)];
params = [ getNormParams(bestParams, bounds) ; ...
           getNormParams(roundBestParams, bounds); ...
           0.9*ones(1, 22); rand(1, 22)];
[score, acc] = vjWrap(XTrain, YTrain, HaarCascade, params);
