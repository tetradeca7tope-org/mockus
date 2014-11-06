% function G = addGaussKernelGradient(X1, X2, Kernels, decomposition, bws)
% 
%   numGroups = numel(decomposition);
%   n1 = size(X1, 1);
%   n2 = size(X2, 1);
%   D = size(X, 2);
%   G = zeros(D, n2, n1);
% 
%   for k = 1:numGroups
%     coords = decomposition{k};
%     K = Kernels(:,:,k);
%     bw = bws(k);
%     X1k = X1(coords, :);
%     X2k = X2(coords, :);
%     G(coords, :,:) = (bw, X1k, X2k, K);
%   end
% end
% 
