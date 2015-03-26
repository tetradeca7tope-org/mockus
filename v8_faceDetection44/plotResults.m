function [ ] = plotResults(rfname)

addpath('results/');

load(rfname);

[numE, numI, numD] = size(boAddCumRewards);
labels = cell(numD,1);
cc = varycolor(numD);

for i=1:numE
  figure(i)
  for j = 1:numD
      plot([1:numI],boAddCumRewards(i,:,j),'Marker','+','color',cc(j,:));
      labels{j} = sprintf('d = %d',numDimsPerGroupCands(j));
      hold on
  end

  titleN = sprintf('Exp:%d    GuessRange:%d   numDiRectEval:%d',i, guessRange,numDiRect);
  title(titleN);
 
  legend(labels,'FontSize',15,'Location','southeast');  
  set(gca,'fontsize',15);
  
  fname = sprintf('figs/%s-%d',rfname,i);
  saveas(i,fname,'png');
end

end
