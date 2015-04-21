% assume the data is already loaded

% compare know decomposition, unknown grouping, and d varies across iterarions

mks = {'+','o','*','x','s','d','^','v','>','<','p','h'};

numExps = size(boUDCumRegrets,1)
totalPts = size(boUDCumRegrets,2);

labels = cell(1+numel(numDimsPerGroupCands),1);

cls = varycolor(1+numel(numDimsPerGroupCands));

boAddSimpleRegrets = mean(boAddSimpleRegrets,1);

for i=1:numel(numDimsPerGroupCands)
  semilogy((1:totalPts),boAddSimpleRegrets(1,:,i),'color',cls(i,:),'Marker',mks{i});
  labels{i} = sprintf('Additive model, d = %d',numDimsPerGroupCands(i));
  hold on
end

boUDSimpleRegrets = mean(boUDSimpleRegrets,1);
semilogy((1:totalPts),boUDSimpleRegrets(1,:),'color',cls(i+1,:),'Marker',mks{i+1});

labels{i+1} = sprintf('Decomposition varies through iterations');

legend(labels, 'FontSize',15,'Location','northeast');
title('Simple Regrets','FontSize',18);
