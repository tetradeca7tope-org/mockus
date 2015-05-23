% assume the data is already loaded
numExps = size(boUDCumRegrets,1)
totalPts = size(boUDCumRegrets,2);
NUM_SHOW = totalPts / 4;
allIters = (1:totalPts);

notRandLines = numel(numDimsPerGroupCands) + numel(choosedM);

cls = varycolor(notRandLines);
clsAdd = cls(1:numel(numDimsPerGroupCands),:);
clsUD = cls(end-numel(choosedM)+1:end,:);

mksAdd = {'^','v','<','>','.'};
mksUD = {'o','p','h'};

labels = cell(notRandLines+1,1);
qqq = round(linspace(1,totalNumQueries, NUM_SHOW));

boAddSimpleRegrets = mean(boAddSimpleRegrets,1);
for i=1:numel(numDimsPerGroupCands)
  semilogy(allIters(qqq), boAddSimpleRegrets(1,qqq,i),'color',clsAdd(i,:),'Marker',mksAdd{i});
  labels{i} = sprintf('Add-d = %d',numDimsPerGroupCands(i));
  hold on
end

boUDSimpleRegrets = mean(boUDSimpleRegrets,1);
for i=1:numel(choosedM)
  semilogy(allIters(qqq),boUDSimpleRegrets(1,qqq,i),'color', clsUD(i,:),'Marker', mksUD{i});
  labels{numel(numDimsPerGroupCands)+i} = choosedM{i};
  hold on
end

randSimpleRegrets = mean(randSimpleRegrets,1);
semilogy(allIters(qqq),randSimpleRegrets(1,qqq),'color','red','Marker','+');
labels{end} = sprintf('Random');

legend(labels, 'FontSize',12,'Location','northeast');
ylabel('Simple Regrets','FontSize',12);
xlabel('Number of iterations','FontSize',12);
titleStr = sprintf('(D, d, M) = (%d, %d, %d)', numDims, trueNumDimsPerGroup, floor(numDims/trueNumDimsPerGroup));
title(titleStr,'FontSize',12);
