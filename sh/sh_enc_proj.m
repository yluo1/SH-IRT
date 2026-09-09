function C = sh_enc_proj(max_odr, theta, phi, is_real, enable_disp)
%Spherical harmonic expansion projection function

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
%codegen('sh_enc_proj', '-o', 'sh/sh_enc_proj_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Plot projection function for higher orders

% theta = pi/2;
% phi = pi/3;
% C_1 = sh_enc_proj(1, theta, phi, false, true);
% C_2 = sh_enc_proj(2, theta, phi, false, true);
% C_3 = sh_enc_proj(3, theta, phi, false, true);

arguments
    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;
    theta (:,1) double = [0];
    phi   (:,1) double = [0];
    is_real (1,1) logical = false;
    enable_disp (1,1) logical = false;
end

C = sh_val(max_odr, theta, phi, is_real)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if enable_disp && coder.target('MATLAB')

    %Mercator plot
    sh_plt(C, 'mercator', is_real);

end