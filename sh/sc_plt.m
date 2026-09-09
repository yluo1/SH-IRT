function h_fig = sc_plt(theta, phi, r, options)
%Plot spherical coordinate scatter plot

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%theta:             [N x 1] Co-latitude (radians), (0, pi)
%phi:               [N x 1] Azimuth (radians), (0, 2 * pi)
%r:                 [N x 1] Radius (meters)

%options:           struct

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h_fig:             Handle to figure

arguments
    theta (:,1) double = [0];
    phi   (:,1) double = [0];
    r (:,1) double = [1];

    options.markersize = 36;
    options.fontsize = 16;
end

assert(numel(theta) == numel(phi) && numel(theta) == numel(r), 'theta, phi, r size mismatch')

[x, y, z] = sph2cart(phi, pi/2 - theta, r);

h_fig = figure; 

scatter3(x, y, z, options.markersize, abs((bsxfun(@rdivide, [x, y, z], r) + 1)/2 ) );

xlabel('X Axis', 'fontsize', options.fontsize);
ylabel('Y Axis', 'fontsize', options.fontsize);
zlabel('Z Axis', 'fontsize', options.fontsize);
set(gca, 'fontsize', options.fontsize - 1);