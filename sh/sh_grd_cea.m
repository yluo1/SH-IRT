function [theta, phi, v, h_fig] = sh_grd_cea(name, options)
% Spatial window directions for loudspeaker measurements

% Reference: Consumer Electronics Association (CEA). Standard Method of Measurement 
% for In-Home Loudspeakers. ANSI/CEA-2034-A (ANSI/CTA-2034-A), Feb. 2015.

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%name:      String, window name {'onaxis', 'listening', 'early-reflections'}

%options:               struct
%options.enable_disp:   If true, plot directions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%theta:         [N x 1]  Co-latitude [0, pi]
%phi:           [N x 1]  Azimuth [0, 2 * pi)

%v:             [N x 3]  Cartesian coordinates, unit norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample usage: Generate CEA-2034-A measurement directions

% [theta_OA, phi_OA, v_OA, h_fig_OA]     = sh_grd_cea('on-axis', 'enable_disp', true);
% [theta_lis, phi_lis, v_lis, h_fig_lis] = sh_grd_cea('listening', 'enable_disp', true);
% [theta_ER, phi_ER, v_ER, h_fig_ER]     = sh_grd_cea('early-reflections', 'enable_disp', true);

% exportgraphics(h_fig_lis, fullfile('figs/figs_grd/', 'grd_cea2034a_lis.png'));
% exportgraphics(h_fig_ER, fullfile('figs/figs_grd/', 'grd_cea2034a_ER.png'));

arguments
    name (1,:) char {mustBeMember(name, {'on-axis', 'listening', 'early-reflections'})} = 'listening';
    options.enable_disp (1,1) logical = false;
end

if strcmp(name, 'on-axis')

    azim_elev_deg = [0, 0];

elseif strcmp(name, 'listening')

    azim_elev_deg = [-30 0; -20 0; -10 0; 0 0; 10 0; 20 0; 30 0; 0 -10; 0 10];

elseif strcmp(name, 'early-reflections')

    floor_azim_elev_deg = [0 -20; 0 -30; 0 -40];
    ceil_azim_elev_deg = [0 40; 0 50; 0 60];
    front_azim_elev_deg = [-30 0; -20 0; -10 0; 0 0; 10 0; 20 0; 30 0];
    side_azim_elev_deg = [-40 0; -50 0; -60 0; -70 0; -80 0; 40 0; 50 0; 60 0; 70 0; 80 0];
    rear_azim_elev_deg = [90:10:270; zeros(1, 19)]';

    azim_elev_deg = [floor_azim_elev_deg; ceil_azim_elev_deg; front_azim_elev_deg; side_azim_elev_deg; rear_azim_elev_deg];

else
    error('Unsupported name');
end

phi     = deg2rad(azim_elev_deg(:,1));
theta   = deg2rad(90 - azim_elev_deg(:,2));

if nargout > 2
    N = numel(theta);
    v = zeros(N, 3);
    [v(:, 1), v(:, 2), v(:, 3)] = sph2cart(phi, pi/2 - theta, ones(N, 1));
end

% Plotting
h_fig = [];
if options.enable_disp && coder.target("MATLAB")
    h_fig = sc_plt(theta, phi, ones(size(theta)), 'title_name', ['CEA-2034-A: ', name], 'filled', true); 
    xlim([-1.1 1.1]);
    ylim([-1.1 1.1]);
    zlim([-1.1 1.1]);
    view(45, 30);
    axis square;
end