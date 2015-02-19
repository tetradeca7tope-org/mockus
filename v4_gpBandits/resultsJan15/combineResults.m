
% baseStr
% clearvars -except baseStr;

fileSearchStr = sprintf('%s/*.mat', baseStr);
files = dir(fileSearchStr);
f1 = sprintf('%s/%s', baseStr, files(1).name);
load(f1);


% Now add all of it
for i = 2:numel(files)
  fi = sprintf('%s/%s', baseStr, files(i).name);
  L = load(fi);

  if size(L.boKDSimpleRegrets, 2) == 910
  
  %   sprintf('%s, %d', files(i).name,
    % Now add them one by one
    boKDSimpleRegrets = [boKDSimpleRegrets; L.boKDSimpleRegrets];
    boKDCumRegrets = [boKDCumRegrets; L.boKDCumRegrets];

    boUDSimpleRegrets = [boUDSimpleRegrets; L.boUDSimpleRegrets];
    boUDCumRegrets = [boUDCumRegrets; L.boUDCumRegrets];

    boAddSimpleRegrets = [boAddSimpleRegrets; L.boAddSimpleRegrets];
    boAddCumRegrets = [boAddCumRegrets; L.boAddCumRegrets];

    randSimpleRegrets = [randSimpleRegrets; L.randSimpleRegrets];
    randCumRegrets = [randCumRegrets; L.randCumRegrets];

  end

end

% Remove all zeros
incIdxs = boKDSimpleRegrets(:, end) ~= 0;
boKDSimpleRegrets = boKDSimpleRegrets(incIdxs, :);
boUDSimpleRegrets = boUDSimpleRegrets(incIdxs, :);
boAddSimpleRegrets = boAddSimpleRegrets(incIdxs, :, :);
randSimpleRegrets = randSimpleRegrets(incIdxs, :);
boKDCumRegrets = boKDCumRegrets(incIdxs, :);
boUDCumRegrets = boUDCumRegrets(incIdxs, :);
boAddCumRegrets = boAddCumRegrets(incIdxs, :, :);
randCumRegrets = randCumRegrets(incIdxs, :);
size(boAddSimpleRegrets),

% Now save the results
saveFileName = sprintf('%s.mat', baseStr);

  save(saveFileName, 'numDims', 'trueNumDimsPerGroup', 'func', ...
    'funcProperties', 'trueMaxVal', 'numDimsPerGroupCands',  ...
    'numIters', 'totalNumQueries', 'numdCands', ...
    'boKDHistories', 'boUDHistories', 'boAddHistories', 'randHistories', ...
    'diRectHist', 'diRectHistory', ...
    'boKDSimpleRegrets', 'boUDSimpleRegrets', 'boAddSimpleRegrets', ...
    'randSimpleRegrets', 'diRectSimpleRegret', ...
    'boKDCumRegrets', 'boUDCumRegrets', 'boAddCumRegrets', 'randCumRegrets', ...
    'diRectCumRegret' );
