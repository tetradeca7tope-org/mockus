% Unit test for vjWrap.m

load VJ/vjdata
addpath VJ/
addpath ../utils/
HaarCascade = GetHaarCasade('VJ/HaarCascades/haarcascade_frontalface_alt.mat');

bounds = [ 0 2; 0 10; 5 15; 10 20; 10 20; 15 25; 15 25; 20 30; 20 30; 30 40; ...
  30 40; 45 55; 50 60; 45 55; 65 75; 60 70; 65 75; 75 85; 80 90; 85 95; 100 110; ...
  100 110];
bestParams = [0.8227 6.9566 9.4985 18.4130 15.3241 21.0106 23.9188 24.5279 ... 
  27.1534 34.5541 39.1073 50.6105 54.6201 50.1697 66.6691 67.6989 69.2288 ...
  79.2491 87.6960 90.2533 104.7492 105.7611];
roundBestParams = round(bestParams);
vjOptions.numXTrain = 100;
vjOptions.doubleParams = false;
% params = [  0.9*ones(1, 22); rand(1, 22)];
% params = [ getNormParams(bestParams, bounds) ; ...
%            getNormParams(roundBestParams, bounds); ...
%            0.5*ones(1,22); 0.9*ones(1, 22); rand(5, 22)];
params = rand(500, 22);
[score, acc] = vjWrap(XTrain, YTrain, HaarCascade, params, bounds, vjOptions);
