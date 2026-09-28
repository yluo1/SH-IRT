function [theta, phi] = sh_grd_caps(N_theta, N_phi, options)
%Generate spherical coordinate grid over uniform co-latitude and azimuth

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%N_theta:       Number of uniform points along co-latitude, excluding poles
%N_phi:         Number of uniform points along azimuth

%options:               Struct
%options.inc_poles:     Logical, if true, include co-latitude poles
%options.axis:          String, rotate coordinates onto axis {'z', 'y', 'x'}
%options.enable_disp:   Logical, if true, plot spherical coordinates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%theta:         [N  x 1]  Co-latitude [0, pi]
%phi:           [N  x 1]  Azimuth [0, 2 * pi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate and plot points

%[theta_z, phi_z] = sh_grd_caps(20, 20, 'enable_disp', true);
%[theta_y, phi_y] = sh_grd_caps(20, 20, 'axis', 'y', 'enable_disp', true);
%[theta_x, phi_x] = sh_grd_caps(20, 20, 'axis', 'x', 'enable_disp', true);

arguments
    N_theta (1,1) double {mustBePositive, mustBeInteger} = 10;
    N_phi   (1,1) double {mustBePositive, mustBeInteger} = 10;

    options.inc_poles (1,1) logical = true;
    options.axis (1,1) char {mustBeMember(options.axis, {'z', 'y', 'x'})} = 'z';
    options.enable_disp (1,1) logical = false;
end

theta_list = linspace(0, pi, N_theta + 2);
theta_list = theta_list(2:end-1);

phi_list = linspace(0, 2 * pi, N_phi + 1);
phi_list = phi_list(1:end-1);

[theta_mat, phi_mat] = meshgrid(theta_list, phi_list);

theta = theta_mat(:);
phi = phi_mat(:);

% Include poles
if options.inc_poles
    theta = [theta; 0; pi];
    phi   = [phi; 0; 0];
end

% Rotation
if ~strcmp(options.axis, 'z')

    if strcmp(options.axis, 'x')
        R = roty(90);
    elseif strcmp(options.axis, 'y')
        R = rotx(90);
    else
        R = eye(3);
    end
    
    [x, y, z] = sph2cart(phi, pi/2 - theta, ones(size(theta)));    
    v = [x, y, z] * R';
    
    [az, elev] = cart2sph(v(:, 1), v(:, 2), v(:, 3));
    phi = az;
    theta = pi/2 - elev;

end

% Plotting
if options.enable_disp
    sc_plt(theta, phi, ones(size(theta)));
end