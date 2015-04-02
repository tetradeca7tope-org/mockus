% assume the data is already loaded

% compare know decomposition, unknown grouping, and d varies across iterarions

figure(1)
% simple regrets
plot((1:numIters),boKDSimpleRegrets(1,:));
hold on
for i=1:
plot((1:numIters),boAddSimpleRegrets(1,:));
hold on
plot((1:numIters),boUDSimpleRegrets(1,:));
hold on



figure(2)
% cumulative regrets

figure(3)
% maximum values

