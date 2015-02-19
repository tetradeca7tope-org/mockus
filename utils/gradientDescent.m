function [fmin, xmin] = gradientDescent(f, gradf, init_pt, params)
% Performs gradient descent on the function f. f is a function handle for the
% function to be minimized. gradf is a function handle for its gradient (or any
% descent direction for the function).
% We use back/forward-tracking line search.

  % Parameters for Back/Forward Tracking
  ALPHA = 0.3; % B&V recommends alpha \in (0.01, 0.3)
  INV_BETA = 1.5; % According to B&V BETA \in (0,1)

  % Determine termination criteria
  if ~exist('params', 'var') params = struct();
  end
  if isfield(params, 'num_iters') max_num_iters = params.num_iters;
  else max_num_iters = 100; % don't do more than 100 steps
  end
  if isfield(params, 'norm_tol') norm_tol = params.norm_tol;
  else norm_tol = 1e-15;
  end
  if isfield(params, 'init_step_size') init_step_size = params.init_step_size;
  else init_step_size = 1;
  end

  % Initialize Gradient Descent
  iter = 0;
  diff = inf;
  curr_pt = init_pt;
  prev_val = inf;
  curr_val = f(init_pt);

  all_vals = zeros(max_num_iters, 1);
  for iter = 1:max_num_iters
    % Prelims
    prev_val = curr_val;

    % Determine whether to do forward or backward tracking
    gradx = gradf(curr_pt);
    norm_gradx = norm(gradx);
    if norm_gradx < norm_tol
      break; % that means you are very close to the minimum
    end
    desc_dir = - init_step_size * gradx / norm(gradx); % scale by init_step_size
    dec_val = gradx' * desc_dir;
    if f(curr_pt + desc_dir) >= curr_val + ALPHA * dec_val
      % If in your first prediction the point on the line is below your function
      % use Backtracking to choose a smaller step size
      step_size = backTracking(f, curr_pt, curr_val, desc_dir, dec_val, ...
                               ALPHA, INV_BETA);
    else
      % Otherwise your initial step size could be too small so use forwrad
      % tracking to choose a larger step size.
      step_size = fwdTracking(f, curr_pt, curr_val, desc_dir, dec_val, ...
                              ALPHA, INV_BETA);
    end

    % Finally perform the updates
    curr_pt = curr_pt + step_size * desc_dir;
    curr_val = f(curr_pt);
    all_vals(iter) = curr_val;
%     curr_pt, curr_val, step_size,
  end % end while

  fmin = curr_val;
  xmin = curr_pt;
end % end of function


function step_size = backTracking(f, curr_pt, curr_val, desc_dir, dec_val, ...
  alpha, inv_beta) 
  step_size = 1/inv_beta;
  while f(curr_pt + step_size*desc_dir) > curr_val + alpha*step_size*dec_val
    step_size = step_size/ inv_beta;
  end
end


function step_size = fwdTracking(f, curr_pt, curr_val, desc_dir, dec_val, ...
  alpha, inv_beta) 
  step_size = inv_beta;
  while f(curr_pt + step_size*desc_dir) < curr_val + alpha*step_size*dec_val
    step_size = step_size * inv_beta;
  end
  step_size = step_size / inv_beta;
end

