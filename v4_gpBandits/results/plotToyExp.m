% assume the data is already loaded

% compare know decomposition, unknown grouping, and d varies across iterarions

totalPts = size(boKDCumRegrets,2);

figure(1)
% simple regrets
labels = cell(2+numel(numDimsPerGroupCands),1);

semilogy((1:totalPts),boKDSimpleRegrets(2,:),'color','green','Marker','+');
labels{1} = sprintf('Known Decomposition: true d = %d', trueNumDimsPerGroup);
hold on

for i=1:numel(numDimsPerGroupCands)
  i=2;
  tmp = numel(numDimsPerGroupCands);
  semilogy((1:totalPts),boAddSimpleRegrets(2,:,i),'color',[1 i/tmp 1-i/tmp], 'linewidth',3);
  labels{i+1} = sprintf('Additive model, d = %d',numDimsPerGroupCands(i));
  hold on
end

semilogy((1:totalPts),boUDSimpleRegrets(2,:),'color','blue','Marker','o');


labels{i+2} = sprintf('Decomposition varies through iterations');

%legend(labels, 'FontSize',15,'Location','northeast');
title('Simple Regrets','FontSize',18);
