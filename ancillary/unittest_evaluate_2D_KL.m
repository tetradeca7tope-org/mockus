% unit test for estimate_2D_KL using 2 gaussians. 

% Define the 2 Gaussians
m0 = [0; 0];
S0 = [1 0; 0 1];
m1 = [2; 1];
S1 = [1 0.1; 0.1 2];

% evaluate analytic KL
ana_kl = 0.5* (trace(S1\S0) + (m1-m0)'*(S1\(m1-m0)) ...
               - 2 - log(det(S0)/det(S1)) );

% Now evaluate numerically
N = 1000;
x = linspace(-4, 6, N);
y = linspace(-8, 10, N);
[X, Y] = meshgrid(x,y);
XX = [X(:), Y(:)];
logP0 = log( mvnpdf(XX, m0', S0) );
logP1 = log( mvnpdf(XX, m1', S1) );
num_kl = estimate_2D_KL(XX, logP0, logP1);

% print out results
fprintf('Analytical KL: %f\nNumerical KL: %f\n', ana_kl, num_kl);
