function [theta, phi, v] = sh_grd_fib(N)
%Generate Fibonacci sphere points (~uniform spaced on unit sphere)

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%N:             Number of points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%theta:         [N x 1]  Co-latitude [0, pi]
%phi:           [N x 1]  Azimuth [0, 2 * pi)

%v:             [N x 3]  Cartesian coordinates, unit norm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_grd_fib', '-o', 'sh/sh_grd_fib_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate and plot points

%N = 1000;
%[theta, phi] = sh_grd_fib(N);
%sc_plt(theta, phi, ones(size(theta)));

arguments
    N (1,1) double {mustBePositive, mustBeInteger} = 100;    
end

golden_ratio = (1 + sqrt(5)) / 2;   

n = (0:(N - 1))';

theta = acos(1 - 2 * (n + 0.5) / N); %Avoids concentration of points at poles
phi = mod(2 * pi * n / golden_ratio, 2 * pi);

if nargout > 2
    v = zeros(N, 3);
    [v(:, 1), v(:, 2), v(:, 3)] = sph2cart(phi, pi/2 - theta, ones(N, 1));
end