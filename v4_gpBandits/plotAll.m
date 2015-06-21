
clear all;
dirName = 'resultsJan15';
dirName = 'resultsFeb26';
dirSearchName = sprintf('%s/*.mat', dirName);
files = dir(dirSearchName);

% Now add all of it
for i = 1:numel(files)
  fi = sprintf('%s/%s', dirName, files(i).name);
%   L = load(fi);
  saveFileName = fi;
  plotGPBResultsDuplicate;
  pause;

end


