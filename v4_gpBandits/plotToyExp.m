% assume the data is already loaded
mks = {'^','v','<','>','p','h','.'};

numExps = size(boUDCumRegrets,1)
totalPts = size(boUDCumRegrets,2);
NUM_SHOW = totalPts / 4;
allIters = (1:totalPts);

labels = cell(2+numel(numDimsPerGroupCands),1);
cls = varycolor(numel(numDimsPerGroupCands));

qqq = round(linspace(1,totalNumQueries, NUM_SHOW));

boAddSimpleRegrets = mean(boAddSimpleRegrets,1);
for i=1:numel(numDimsPerGroupCands)
  semilogy(allIters(qqq), boAddSimpleRegrets(1,qqq,i),'color',cls(i,:),'Marker',mks{i});
  labels{i} = sprintf('Add-d = %d',numDimsPerGroupCands(i));
  hold on
end

boUDSimpleRegrets = mean(boUDSimpleRegrets,1);
semilogy(allIters(qqq),boUDSimpleRegrets(1,qqq),'color','black','Marker','o');
labels{1+i} = sprintf('CHOOSEdM');
hold on

randSimpleRegrets = mean(randSimpleRegrets,1);
semilogy(allIters(qqq),randSimpleRegrets(1,qqq),'color','red','Marker','+');
labels{2+i} = sprintf('Random');

legend(labels, 'FontSize',12,'Location','northeast');
ylabel('Simple Regrets','FontSize',12);
xlabel('Number of iterations','FontSize',12);
titleStr = sprintf('(D, d, M) = (%d, %d, %d)', numDims, trueNumDimsPerGroup, floor(numDims/trueNumDimsPerGroup));
title(titleStr,'FontSize',12);
