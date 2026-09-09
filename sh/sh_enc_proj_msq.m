function C = sh_enc_proj_msq(max_odr, theta, phi, is_real, enable_disp)
%Spherical harmonic expansion magnitude squared projection function

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%max_odr:       Max SH order

%theta:         [N x 1]  Co-latitude (0, pi)
%phi:           [N x 1]  Azimuth (0, 2 * pi)

%is_real:       Logical, if true, evaluate real SH
%enable_disp:   Logical, if true, plot function expansion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:             [(max_odr + 1)^2 x N] SH coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_enc_proj_msq', '-o', 'sh/sh_enc_proj_msq_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Plot magnitude squared projection function for higher orders.
%Compare with direct projections

%is_real = false;
%dB_lim = [-120, 0];

% theta = pi/2;
% phi = pi/3;
% C_proj_msq_2 = sh_enc_proj_msq(2, theta, phi, is_real);
% C_proj_2 = sh_enc_proj(2, theta, phi, is_real);

% C_proj_msq_4 = sh_enc_proj_msq(4, theta, phi, is_real);
% C_proj_4 = sh_enc_proj(4, theta, phi, is_real);

% C_proj_msq_6 = sh_enc_proj_msq(6, theta, phi, is_real);
% C_proj_6 = sh_enc_proj(6, theta, phi, is_real);

%sh_plt(C_proj_msq_2, 'mercator', is_real, 'title_name', 'Magnitude Square Projection, Max-Odr. 2', 'dB_lim', dB_lim);
%sh_plt(C_proj_2, 'mercator', is_real, 'title_name', 'Projection, Max-Odr. 2', 'dB_lim', dB_lim);
%sh_plt(C_proj_msq_4, 'mercator', is_real, 'title_name', 'Magnitude Square Projection, Max-Odr. 4', 'dB_lim', dB_lim);
%sh_plt(C_proj_4, 'mercator', is_real, 'title_name', 'Projection, Max-Odr. 4', 'dB_lim', dB_lim);
%sh_plt(C_proj_msq_6, 'mercator', is_real, 'title_name', 'Magnitude Square Projection, Max-Odr. 6', 'dB_lim', dB_lim);
%sh_plt(C_proj_6, 'mercator', is_real, 'title_name', 'Projection, Max-Odr. 6', 'dB_lim', dB_lim);

arguments
    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;
    theta (:,1) double = [0];
    phi   (:,1) double = [0];
    is_real (1,1) logical = false;
    enable_disp (1,1) logical = false;
end

N = numel(theta);

half_odr = floor(max_odr / 2);
C_proj = sh_val(half_odr, theta, phi, is_real)'; %[floor(max_odr / 2) x N]

C = [sh_msq(C_proj, is_real); zeros((max_odr + 1)^2 - (2 * half_odr + 1)^2, N)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if enable_disp && coder.target('MATLAB')

    %Mercator plot
    sh_plt(C, 'mercator', is_real);

end