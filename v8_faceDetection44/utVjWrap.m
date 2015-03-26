% Unit test for vjWrap.m

load VJ/vjdata
addpath VJ/
addpath ../utils/
HaarCascade = GetHaarCasade('VJ/HaarCascades/haarcascade_frontalface_alt.mat');

bounds = [ 0 2; 0 10; 5 15; 10 20; 10 20; 15 25; 15 25; 20 30; 20 30; 30 40; ...
  30 40; 45 55; 50 60; 45 55; 65 75; 60 70; 65 75; 75 85; 80 90; 85 95; 100 110; ...
  100 110];

upperBounds = bounds .* 2;
bounds = [bounds; upperBounds];

vjOptions.doubleParams = true;

params = rand(1,44);
[score, acc] = vjWrap(XTrain, YTrain, HaarCascade, params, bounds, vjOptions);
