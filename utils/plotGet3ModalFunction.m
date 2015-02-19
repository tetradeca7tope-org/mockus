
clear all;
close all;

[f, fProps] = get3ModalFunctionDuplicate(2);
N = 40;

t = linspace(-1,1,N);
[X,Y] = meshgrid(t,t);
XY = [X(:), Y(:)];

F = f(XY);
Z = reshape(F, N, N);

% mesh(X, Y, Z);
surf(X, Y, Z);
% contour(X, Y,Z);
set(gca, 'xtick', []);
set(gca, 'ytick', []);
