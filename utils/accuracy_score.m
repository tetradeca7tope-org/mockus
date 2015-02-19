function score = accuracy_score(true_labels, cluster_labels)
%ACCURACY Compute clustering accuracy using the true and cluster labels and
%   return the value in 'score'.
%
%   Input  : true_labels    : N-by-1 vector containing true labels
%            cluster_labels : N-by-1 vector containing cluster labels
%
%   Output : score          : clustering accuracy

% Compute the confusion matrix 'cmat', where
%   col index is for true label (CAT),
%   row index is for cluster label (CLS).

  % First of all convert this back to the required format.
  if size(true_labels, 2) > 1
    T = zeros(size(true_labels,1), 1);
    for i = 1:size(true_labels,1)
      T(i) = find(true_labels(i,:));
    end
    true_labels = T;
  end
  % Do the same for predicts
  if size(cluster_labels, 2) > 1
    P = zeros(size(cluster_labels,1), 1);
    for i = 1:size(cluster_labels,1)
      P(i) = find(cluster_labels(i,:));
    end
    cluster_labels = P;
  end

% Print out the number of instances detected per class labels.
unique_labels = unique(cluster_labels);
for i = 1:numel(unique_labels)
  fprintf('Class lbl #%d : %d\n', unique_labels(i), ...
          sum(cluster_labels == unique(i)));
end

n = length(true_labels);
cat = spconvert([(1:n)' true_labels ones(n,1)]);
cls = spconvert([(1:n)' cluster_labels ones(n,1)]);
cls = cls';
cmat = full(cls * cat);

%
% Calculate accuracy
%
[match, cost] = hungarian(-cmat);
score = 100*(-cost/n);
