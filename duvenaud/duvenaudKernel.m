function KXY = duvenaudKernel(X, Y, baseKernel, maxOrder, scales)

  numDims = size(X, 2);
  numX = size(X, 1);
  numY = size(Y, 1);
  
  if isfloat(baseKernel)
    % Then baseKernel should be the bandwidth of a Gauss Kernel
    baseKernel = @(x,y) GaussKernel(baseKernel, x, y);
  end

  if numel(scales) == 1
    scales = scales * ones(maxOrder, 1);
  end
  
  % We need the power Sums
  dimBaseKs = zeros(numX, numY, numDims);
  for i = 1:numDims
    dimBaseKs(:,:,i) = baseKernel( X(:,i), Y(:,i) );
  end
  powerSums = zeros(numX, numY, maxOrder)
  for a = 1:maxOrder
    powerSums(:,:,i) = sum(dimBaseKs.^a , 3); ;
  end

  % Now we construct the additive Kernels of each order
  orderAddKernels = zeros(numX, numY, maxOrder + 1);
  orderAddKernels(:,:,1) = ones(numX, numY);
  for d = 1:maxOrder
    accum = zeros(numX, numY);
    for a = 1:d
      accum = accum + (-1)^(a-1)*orderAddKernels(:,:,d+1-a).*powerSums(:,:,a);
    end
    orderAddKernels(:,:,d+1) = 1/d * accum; 
  end

  % Now compute the kernel
  KXY = zeros(numX, numY);
  for a = 1:n
    KXY = scales(a) * orderAddKernels(:,:,a+1);
  end

end
