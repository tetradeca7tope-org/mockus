function plot1DFunction(f, range, color, desc)

  t = linspace(range(1), range(2), 200)';
  ft = f(t);

  plot(t, ft, desc, 'Color', color);

end
