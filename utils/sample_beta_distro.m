% params
PLOT_OK = 0;

% Plot the function first
if PLOT_OK
  th = linspace(0,1);
  f1 = figure; plot(th, th.^(alp-1) .* (1 - th).^(bep-1) / beta(alp, bep) );
  title('prior');
end

% Sample p from this prior
theta = dirichlet_sample([alp, bep]); p = theta(1);
fprintf('p = %f\n', p);
% Now generate N Bernoulli points from this p
X = double(rand(N,1) < p);
sumX = sum(X);
fprintf('Number of positive trials : %d\n', sumX);

