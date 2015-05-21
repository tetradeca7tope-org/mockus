% assume the data is already loaded

mks = {'+','^','v','<','>','p','h','.'};

numExps = size(boUDCumRegrets,1)
totalPts = size(boUDCumRegrets,2);

labels = cell(2+numel(numDimsPerGroupCands),1);
cls = varycolor(2+numel(numDimsPerGroupCands));

boAddSimpleRegrets = mean(boAddSimpleRegrets,1);

for i=1:numel(numDimsPerGroupCands)
  semilogy((1:totalPts),boAddSimpleRegrets(1,:,i),'color',cls(i,:),'Marker',mks{i});
  labels{i} = sprintf('Add-d = %d',numDimsPerGroupCands(i));
  hold on
end

boUDSimpleRegrets = mean(boUDSimpleRegrets,1);
semilogy((1:totalPts),boUDSimpleRegrets(1,:),'color',cls(1+i,:),'Marker',mks{1+i});
labels{1+i} = sprintf('Decomposition varies');
hold on

randSimpleRegrets = mean(randSimpleRegrets,1);
semilogy((1:totalPts),randSimpleRegrets(1,:),'color',cls(2+i,:),'Marker',mks{2+i});
labels{2+i} = sprintf('Random optimization');

legend(labels, 'FontSize',15,'Location','northeast');
title('Simple Regrets','FontSize',18);
