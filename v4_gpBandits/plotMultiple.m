for i=1:9
  mks = {'+','o','*','x','s','d','^','v','>','<','p','h'};
  fileStr = sprintf('results/toPlot/%d',i);
  load(fileStr);

  subplot(3,3,i);
  numExps = size(boUDSimpleRegrets,1);
  totalPts = size(boUDSimpleRegrets,2);

  labels = cell(1+numel(numDimsPerGroupCands),1);
  cls = varycolor(1+numel(numDimsPerGroupCands));

  boAddSimpleRegrets = mean(boAddSimpleRegrets,1);
  for i=1:numel(numDimsPerGroupCands)
    semilogy((1:totalPts),boAddSimpleRegrets(1,:,i),'color',cls(i,:),'Marker',mks{i});
    labels{i} = sprintf('Add-d=%d',numDimsPerGroupCands(i));
    hold on
  end

  boUDSimpleRegrets = mean(boUDSimpleRegrets,1);
  semilogy((1:totalPts),boUDSimpleRegrets(1,:),'color','black','Marker',mks{i+1});
  labels{i+1} = sprintf('Decomposition varies');

  legend(labels,'FontSize',12,'Location','northeast');
  ylabel('Simple Regrets','FontSize',12);
  xlabel('Number of iterations','FontSize',12);
  titleStr = sprintf('(D, d, M) = (%d, %d, %d)', numDims, trueNumDimsPerGroup, floor(numDims/trueNumDimsPerGroup));
  title(titleStr,'FontSize',12);
end
