function plot2DFunction(f, range, plotType, plotArgs)
% f is a function handle taking in an a nx2 array and returning an nx1 output.
% range specifies the 2D range. If empty taken to be [0 1 0 1]
% plotType is either 'mesh' or 'contour'

  if isempty(range)
    range = [0 1 0 1];
  end
  if isequal(size(range), [2 2])
    range = [range(1,:), range(2,:)];
  end

  if ~exist('plotArgs', 'var') | isempty(plotArgs)
    plotArgs = '';
  end

  res = 100;
  t1 = linspace(range(1), range(2), res)';
  t2 = linspace(range(3), range(4), res)';
  [T1, T2] = meshgrid(t1, t2);
  th = [T1(:) T2(:)];
  fth = f(th);
  FTH = reshape(fth, res, res);
  if strcmp(plotType, 'mesh')
    mesh(T1, T2, FTH);
  elseif strcmp(plotType, 'contour')
    contour(T1, T2, FTH, plotArgs);
  elseif strcmp(plotType, 'contour3')
    contour3(T1, T2, FTH, 'contour');
  else
    error('Unknown plot type');
  end
  
end
