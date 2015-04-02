% assume the data is already loaded

% compare know decomposition, unknown grouping, and d varies across iterarions

totalPts = size(boKDCumRegrets,2);

figure(1)
% simple regrets
labels = cell(2+numel(numDimsPerGroupCands),1);

plot((1:totalPts),boKDSimpleRegrets(1,:),'color','green','Marker','+');
labels{1} = sprintf('Known Decomposition: true d = %d', trueNumDimsPerGroup);
hold on

for i=1:numel(numDimsPerGroupCands)
  tmp = numel(numDimsPerGroupCands);
  plot((1:totalPts),boAddSimpleRegrets(1,:,i),'color',[1 i/tmp 1-i/tmp], 'linewidth',3);
  labels{i+1} = sprintf('Additive model, d = %d',numDimsPerGroupCands(i));
  hold on
end

plot((1:totalPts),boUDSimpleRegrets(1,:),'color','blue','Marker','o');
hold on
labels{i+2} = sprintf('Decomposition varies through iterations');

legend(labels, 'FontSize',15,'Location','northeast');
title('Simple Regrets','FontSize',18);
saveas(1,'simple');

figure(2)
% cumulative regrets
plot((1:totalPts),boKDCumRegrets(1,:),'color','green','Marker','+');
hold on

for i=1:numel(numDimsPerGroupCands)
  tmp = numel(numDimsPerGroupCands);
  plot((1:totalPts),boAddCumRegrets(1,:,i),'color',[1 i/tmp 1-i/tmp], 'linewidth',3);
  labels{i+1} = sprintf('Additive model, d = %d',numDimsPerGroupCands(i));
  hold on
end

plot((1:totalPts),boUDCumRegrets(1,:),'color','blue','Marker','o');
hold on
legend(labels, 'FontSize',15,'Location','northeast');
title('Cumulative Regrets','FontSize',18);
saveas(2,'cum');

figure(3)
% maximum values
plot((1:totalPts),boKDHistories(1,:),'color','green','Marker','+');
hold on

for i=1:numel(numDimsPerGroupCands)
  tmp = numel(numDimsPerGroupCands);
  plot((1:totalPts),boAddHistories(1,:,i),'color',[1 i/tmp 1-i/tmp], 'linewidth',3);
  labels{i+1} = sprintf('Additive model, d = %d',numDimsPerGroupCands(i));
  hold on
end

plot((1:totalPts),boUDHistories(1,:),'color','blue','Marker','o');
hold on
legend(labels, 'FontSize',15,'Location','southeast');
title('Maximum values','FontSize',18);
saveas(3,'maxVals');

