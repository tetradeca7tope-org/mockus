% assume the data is already loaded

% compare know decomposition, unknown grouping, and d varies across iterarions

mks = {'+','o','*','x','s','d','^','v','>','<','p','h'};

totalPts = size(boKDCumRegrets,2);
labels = cell(2+numel(numDimsPerGroupCands),1);

cls = varycolor(2+numel(numDimsPerGroupCands));

semilogy((1:totalPts),boKDSimpleRegrets(1,:),'color',cls(1,:),'Marker',mks{1});
labels{1} = sprintf('Known Decomposition: true d = %d', trueNumDimsPerGroup);
hold on

for i=1:numel(numDimsPerGroupCands)
  tmp = numel(numDimsPerGroupCands);
  semilogy((1:totalPts),boAddSimpleRegrets(1,:,i),'color',cls(i+1,:),'Marker',mks{i+1});
  labels{i+1} = sprintf('Additive model, d = %d',numDimsPerGroupCands(i));
  hold on
end

semilogy((1:totalPts),boUDSimpleRegrets(1,:),'color',cls(i+2,:),'Marker',mks{i+2});

labels{i+2} = sprintf('Decomposition varies through iterations');

legend(labels, 'FontSize',15,'Location','northeast');
title('Simple Regrets','FontSize',18);
