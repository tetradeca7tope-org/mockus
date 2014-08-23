function A = numerical_1D_integration(z, f)

if size(z,2) ~= size(f,2)
  if size(z,2) == 1,
    z = repmat(z,1,size(f,2));
  elseif size(f,2) == 1
    f = repmat(f,1,size(z,2));
  else
    error('sizes of f & z don"t match');
  end
end

m = size(z,1);
n = size(z,2);
A = zeros(1, n);

for i = 1:(m-1)
  accum = (f(i,:) + f(i+1,:)) / 2 .* (z(i+1,:) - z(i,:));
  if abs(accum) < inf
    A = A + accum;
  end
end

end
