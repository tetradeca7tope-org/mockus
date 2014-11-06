% ut for Gauss Kernel Gradient.
% The gradient is partial K/ partial x_i = -K * (z_i - z_j)/h^2

X = 0.1 * [1:4];
Y = [X; 0.1 *[5:8; 9:2:15]];
h = 0.5;

K = GaussKernel(h, X, Y),

G1 = GaussKernelGradient(h, X, Y),
G2 = GaussKernelGradientCoords(h, X, Y, K, []),
G2 = GaussKernelGradientCoords(h, Y, Y, GaussKernel(h, Y, Y), []),

