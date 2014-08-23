function [ avg_diagonal, conf_mat_frac, conf_mat, ...
           actual_labels, predicted_labels, ...
           conf_mat_frac_print, conf_mat_print] = ...
            confusion_matrix( Y, predicts, actual_labels)
  % This function computes the confusion matrix for a given classification.
  % Y is a vector containing the actual labels and predicts the predicted
  % labels.
  % This function returns the following:
  % avg_diagonal: The average of the diagonal of the confusion matrix
  % conf_mat_frac: The confusion matrix where the (i,j)th element
  %                represents the fraction of instances carrying the (i)th
  %                label classified to the (j)th label.
  % conf_mat: The confusion matrix as explained above but with the actual
  %           number of elements instead of the fraction
  % actual_labels: A list of the actual labels in Y in the order in which
  %                they appear along the rows of conf_mat.
  % predicted_labels: A list of the predicted labels in predicts in the
  %                   order they appear along the columns of conf_mat.
  % conf_mat_frac_print, conf_mat_frac: Printable forms of conf_mat and
  %                                     conf_mat_frac.
  
  % First of all convert this back to the required format.
  if size(Y, 2) > 1
    YY = zeros(size(Y,1), 1);
    for i = 1:size(Y,1)
      YY(i) = find(Y(i,:));
    end
    Y = YY;
  end
  % Do the same for predicts
  if size(predicts, 2) > 2
    PP = zeros(size(predicts,1), 1);
    for i = 1:size(predicts,1)
      PP(i) = find(predicts(i,:));
    end
    predicts = PP;
  end

%   [Y, predicts],
  
  % These vectors store the unique labels.
  if nargin < 3
    actual_labels = unique(Y);
    predicted_labels = unique(predicts);
  else
    actual_labels = unique(actual_labels);
    predicted_labels = unique(actual_labels);
  end
  num_actual_labels = size(actual_labels, 1);
  num_predicted_labels = size(predicted_labels, 1);
  num_min_labels = min(num_actual_labels, num_predicted_labels);
  
  % Construct the confusion matrix
  conf_mat = zeros(num_actual_labels, num_predicted_labels);
  for a_l = 1 : num_actual_labels
    for p_l = 1 : num_predicted_labels
      union = (Y == actual_labels(a_l)) .* ...
              (predicts == predicted_labels(p_l));
      conf_mat(a_l, p_l) = sum(union);
    end
  end
  
  num_insts_on_each_a_l = sum(conf_mat, 2);
  best_grouping_fracs = max(conf_mat, [], 2) ./ num_insts_on_each_a_l;
  [ ~, sorted_frac_indices] = sort(best_grouping_fracs);
  sorted_actual_labels = actual_labels(sorted_frac_indices);
  
  % Now reorder the confusion matrix so that the maximum values occur along
  % the diagonals
  for i = 1:num_min_labels
    label_ind = find(actual_labels == sorted_actual_labels(i));
    a_l_grouping = conf_mat(label_ind,:);
    [~, match_pl_index] = max(a_l_grouping);
    
    % Swap the 2 columns in the confusion matrix
    temp = conf_mat(:, match_pl_index);
    conf_mat(:, match_pl_index) = conf_mat(:, label_ind);
    conf_mat(:, label_ind) = temp;
    
    % Swap the 2 predicted labels in the confusion matrix
    temp = predicted_labels(match_pl_index);
    predicted_labels(match_pl_index) = predicted_labels(label_ind);
    predicted_labels(label_ind) = temp;
  end
  
  % Now compute the confusion matrix so that the elements contain the
  % fractions of the total assigned to the actual label.
  num_insts_on_each_actual_label = ...
    repmat(num_insts_on_each_a_l, 1, num_predicted_labels);
  conf_mat_frac = conf_mat ./ num_insts_on_each_actual_label;
  
  % Compute the average of the diagonal of the confusion matrix
  avg_diagonal = mean(diag(conf_mat_frac));
  
  conf_mat_print = PadLabels(conf_mat, actual_labels, predicted_labels);
  conf_mat_frac_print = PadLabels(conf_mat_frac, actual_labels, ...
                                  predicted_labels);
  
end


function mat_print = PadLabels(conf_mat, a_labels, p_labels)
% This matrix returns conf_mat with the label indices marked.
  tab_ch = inf;
  mat_print = [p_labels'; tab_ch*ones(1, size(p_labels,1)); conf_mat];
  a_label_col = [tab_ch; tab_ch; a_labels];
  col_2 = tab_ch*ones(size(a_labels,1) + 2, 1);
  mat_print = [a_label_col, col_2, mat_print];
end
