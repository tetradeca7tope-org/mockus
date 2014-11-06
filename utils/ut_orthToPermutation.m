% Unit test for getHungarainCostMatrixOrth.m

clear all;
close all;
clc; clc;

addpath ~/libs/kky-matlab/utils/
addpath ~/libs/kky-matlab/ancillary/

pwr = 1;
D = 1000;   
p = round(D/1);
  
A1 = [ 0.9 0.1 0; 0.02 0.02 1; 0 1 0];

% Now a random matrix
b = randperm(p);
A = zeros(D, p);
for i = 1:p
  A(b(i), i) = 1;
end
A = A + 0.6*randn(D,p);
[A, ~] = qr(A, 0);

for costIdx = 1:2

  % First the matrix A
  C1 = getHungarianCostMatrixOrth(A1, costIdx, pwr);
  P1 = orthToPermutation(A1, costIdx, pwr);
  err1 = mean(P1 ~= [1 3 2]);

  tic,
  C = getHungarianCostMatrixOrth(A, costIdx, pwr);
  toc,
  tic,
  P = orthToPermutation(A, costIdx, pwr);
  toc,
  err = mean(P ~= b);

  fprintf('costIdx %d: err1: %0.5f err2: %0.5f\n', costIdx, err1, err);

end

