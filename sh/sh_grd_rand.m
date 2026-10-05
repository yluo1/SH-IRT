function [theta, phi, v] = sh_grd_rand(N)
%Sample uniform directions over unit sphere

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%N:             Number of points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%theta:         [N x 1]  Co-latitude [0, pi]
%phi:           [N x 1]  Azimuth [0, 2 * pi)

%v:             [N x 3]  Cartesian coordinates, unit norm

theta = acos(2 * rand(N, 1) - 1);
phi = rand(N, 1) * 2 * pi;

if nargout > 2
    v = zeros(N, 3);
    [v(:, 1), v(:, 2), v(:, 3)] = sph2cart(phi, pi/2 - theta, ones(N, 1));
end