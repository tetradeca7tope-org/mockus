% assume the data is already loaded

% compare know decomposition, unknown grouping, and d varies across iterarions

mks = {'+','o','*','x','s','d','^','v','>','<','p','h'};

numExps = size(boUDCumRegrets,1)
totalPts = size(boUDCumRegrets,2);

%labels = cell(2+numel(numDimsPerGroupCands),1);
lables = cell(2);

cls = varycolor(2+numel(numDimsPerGroupCands));
cls = varycolor(2);

% boAddSimpleRegrets = mean(boAddSimpleRegrets,1);
% 
% for i=1:numel(numDimsPerGroupCands)
%   semilogy((1:totalPts),boAddSimpleRegrets(1,:,i),'color',cls(i,:),'Marker',mks{i});
%   labels{i} = sprintf('Add-d = %d',numDimsPerGroupCands(i));
%   hold on
% end

boUDSimpleRegrets = mean(boUDSimpleRegrets,1);
semilogy((1:totalPts),boUDSimpleRegrets(1,:),'color',cls(1,:),'Marker',mks{1});
labels{1} = sprintf('Decomposition varies');
hold on

randSimpleRegrets = mean(randSimpleRegrets,1);
semilogy((1:totalPts),randSimpleRegrets(1,:),'color',cls(2,:),'Marker',mks{2});
labels{2} = sprintf('Random optimization');


legend(labels, 'FontSize',15,'Location','northeast');
title('Simple Regrets','FontSize',18);
