function [x,y]=OneScaleObjectDetection( x, y, Scale, IntegralImages, w,h, ...
HaarCasade, thresholds, ithImage)

% [x,y]=OneScaleObjectDetection( x, y, Scale, IntegralImages, w,h,HaarCasade)
%

% Calculate the mean 
InverseArea = 1 / (w*h);
mean =  GetSumRect(IntegralImages.ii,x,y,w,h)*InverseArea;
% Use the mean and squared integral image to calculate the grey-level
% Varianceiance, of every search window
Variance = GetSumRect(IntegralImages.ii2,x,y,w,h)*InverseArea - (mean.^2);
% Convert the Varianceiation to Standard Deviation
Variance(Variance<1)=1; StandardDeviation =sqrt(Variance);

% The haarcasade contains a row of classifier-trees. The classifiers
% are executed one at the time. If a coordinate doesn't pass the classifier
% threshold it is removed, otherwise it goes into the next classifier

% Loop through all classifier stages
% length(HaarCasade.stages),

numStages = length(HaarCasade.stages);
assert(numStages*2 == size(thresholds,1), 'stages and threshods does not match');

info = zeros(1,numStages*2)-1;

for i_stage = 1:numStages,
    stage = HaarCasade.stages(i_stage);
    Trees=stage.trees;
    StageSum = zeros(size(x));

    % Loop through a classifier tree
    for i_tree=1:length(Trees)
        Tree = Trees(i_tree).value;
        % Executed the classifier
        TreeSum=TreeObjectDetection(zeros(size(x)),Tree,Scale,x,y,IntegralImages.ii,StandardDeviation,InverseArea);
        StageSum = StageSum + TreeSum;
    end
    % If the StageSum of a coordinate is lower than the treshold it
    % is removed, otherwise it goes into the next stage
    
    threshold = thresholds(i_stage);  
    upperThreshold = thresholds(i_stage + numStages);
    
    %     threshold = stage.stage_threshold; 
    check=StageSum < threshold;  
    %     StageSum, threshold, check,

    % Remove coordinates which don't contain an object 
    x=x(~check);  
    % All coordinates failed
    % on this Scale to detect an
    % object in the image

    if(isempty(x)) 
      % fprintf('Early Stop\n');
      break; 
    end 
    
    y=y(~check);
    StandardDeviation=StandardDeviation(~check);

    if (max(StageSum) > upperThreshold)
%      fprintf('Early Decide\n');
      break;
    end
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    featureInfoSum = sum(StageSum);
%    featureInfoMax  = max(StageSum);
%
%    info(i_stage) = featureInfoSum;
%    info(i_stage + numStages) = featureInfoMax;   
 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('allInfo','allInfo');
% allInfo(ithImage,:) = info;
% save('allInfo','allInfo');

end
