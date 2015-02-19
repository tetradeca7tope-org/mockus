function [labels] = spectralCluster(A, num_clusters)
% Performs Spectral Clustering on the similarity matrix A.

  n = size(A, 1);

  % Construct Normalized Graph Laplacian
  Dsqrt = diag(sqrt(sum(A)));
  Ln = eye(n) - Dsqrt * A * Dsqrt;

  % Obtain the first num_clusters eigenvectors
  % doing full eigendecomposition.
  [V, D] = eig(Ln);
  eigvals = diag(D);
  [~, sort_idxs] = sort(eigvals, 'ascend');
  sort_idxs = sort_idxs(1:num_clusters);
  V = V(:, sort_idxs); D = D(:, sort_idxs);
  Z = bsxfun(@rdivide, V, sqrt(sum(V.^2, 2)) );

  % Run K-means
  init_centres = njw01Init(Z, num_clusters);
  labels = kmeans(Z, [], 'Start', init_centres);

end

function initial_centres = njw01Init(data, num_clusters)
% We need to initialize the centres for K-means clustering.
% For this we adopt the strategy outlined in "On Spectral Clustering :
% Analysis and an Algorithm", Ng, Jordan, Weiss 2001. page 6/8.

  m = size(data, 1); % Num data

  % Define the following vector to hold the centres.
  initial_centres = zeros(num_clusters, num_clusters);

  % Define a num_clusters-1 x m matrix to store the cosine similarities.  
  cosine_similarities = zeros(num_clusters-1, m);

  % First define the following vector holding the L2 norm of all data points.
  % This vector is repeatedly used to compute the cosine similarities.
  data_n = sqrt( sum(data.^2, 2) );

  % Set the first centre to a randomly selected point.
  initial_centres(1,:) = data(randi(m), :);

  % In this loop, on each iteration we shall compute the cosine similarities of
  % each data point to the previous centre.
  % Then we pick the point which has the lowest minimum cosine similarity.
  for i = 2:num_clusters

    cosine_similarities(i-1, :) = abs( (data * initial_centres(i-1, :)' ...
                                        ./ data_n ...
                                        /norm(initial_centres(i-1, :), 2) )' );
    % Now obtain the minimum distances.
    min_dists = max(cosine_similarities(1:i-1, :), [], 1);
    [temp, next_centre_idx] = min(min_dists);
    initial_centres(i, :) = data(next_centre_idx, :);

  end

end

