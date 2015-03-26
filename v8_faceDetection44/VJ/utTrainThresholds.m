% Unit test for trainThresholds

load('vjdata.mat');

addpath SubFunctions/
HaarCascade = GetHaarCasade('HaarCascades/haarcascade_frontalface_alt.mat');
thresholds = 90*ones(22,1);

tic,
score = trainThresholds(XTrain, YTrain, HaarCascade, thresholds);
toc,
