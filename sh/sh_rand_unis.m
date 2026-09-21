function [theta, phi] = sh_rand_unis(N)
%Sample uniform directions over unit sphere

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%N:             Number of points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%theta:         [N x 1]  Co-latitude [0, pi]
%phi:           [N x 1]  Azimuth [0, 2 * pi)

theta = acos(2 * rand(N, 1) - 1);
phi = rand(N, 1) * 2 * pi;
