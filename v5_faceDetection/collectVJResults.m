% Combines the results in all Mat files

dirName = 'resultsComa';
dirSearchStr = sprintf('%s/*.mat', dirName);
files = dir(dirSearchStr);
f1 = sprintf('%s/%s', dirName, files(1).name);
load(f1)

% Now add all of it
for i = 2:numel(files)
  fi = sprintf('%s/%s', dirName, files(i).name);
  L = load(fi);

  % Now add them one by one
  boAddMaxVals = [boAddMaxVals; L.boAddMaxVals];
  randMaxVals = [randMaxVals; L.randMaxVals];

end

% Now save the results
saveFileName = sprintf('%s/combinedResults.mat');
save(saveFileName, 'numDimsPerGroupCands', 'numIters', 'numDims', 'numdCands', ...
  'boAddMaxVals', 'randMaxVals', 'diRectMaxVals', 'totalNumQueries' ...
  );

