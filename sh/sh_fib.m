function [theta, phi] = sh_fib(N)
%Generate Fibonacci sphere points (~uniform spaced on sphere)

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%N:             Number of points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%theta:         [N x 1]  Co-latitude [0, pi]
%phi:           [N x 1]  Azimuth [0, 2 * pi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_fib', '-o', 'sh/sh_fib_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate and plot points

%N = 1000;
%[theta, phi] = sh_fib(N);
%sc_plt(theta, phi, ones(size(theta)));

arguments
    N (1,1) double {mustBePositive, mustBeInteger} = 100;    
end

golden_ratio = (1 + sqrt(5)) / 2;   

n = (0:(N - 1))';

theta = acos(1 - 2 * (n + 0.5) / N); %Avoids concentration of points at poles
phi = mod(2 * pi * n / golden_ratio, 2 * pi);