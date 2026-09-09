function D = sh_rot(C, theta_zyx, is_real, rot_intrinsic)
% Rotate spherical harmonic expansion coefficients by zyx (yaw-pitch-roll)
% Tait-Bryan rotation, intrinsic

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:             [(P+1)^2 x M] SH coefficients (max order P of M number of functions)
%theta_zyx:     [1 x 3]  Radians, counter clock-wise rotation about [+x, +y, +z] axis
%is_real:       Logical, if true, evaluate real SH
%rot_intrinsic: Logical, if true, rotation is intrinsic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%D:             [(P+1)^2 x M] SH coefficients, rotated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_rot', '-o', 'sh/sh_rot_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Rotate directional filtered random field by 90 degrees w.r.t. z, y, x axes,
%and compare with mex

%rng(123);
%max_odr = 8;
%is_real = false;
%C = sh_mul(sh_enc_rbf('SqExp', max_odr, pi/2, 0, 1/2, is_real), sh_rand(max_odr, 1, is_real));
%D_z90 = sh_rot(C, deg2rad([90, 0, 0]), is_real);
%D_y90 = sh_rot(C, deg2rad([0, 90, 0]), is_real);
%D_x90 = sh_rot(C, deg2rad([0, 0, 90]), is_real);
%D_z90_y45 = sh_rot(C, deg2rad([90, 45, 0]), is_real);

%dB_lim = [-60, 12];
%sh_plt(C, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Original');
%sh_plt(D_z90, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Rotate 90 Degree Z-axis');
%sh_plt(D_y90, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Rotate 90 Degree Y-axis');
%sh_plt(D_x90, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Rotate 90 Degree X-axis');
%sh_plt(D_z90_y45, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Rotate 90 Degree Z-axis, 45 Degree Y-axis');

%D_z90_mex = sh_rot_mex(C, deg2rad([90, 0, 0]), is_real);
%err_mex = norm(D_z90 - D_z90_mex)

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
    theta_zyx (1,3) double = [0 0 0];
    is_real (1,1) logical = false;  
    rot_intrinsic (1,1) logical = true;
end

[P, M] = size(C);
P = sqrt(P) - 1;
assert(P - floor(P) == 0, 'Invalid size C');

R_list = sh_rot_mat(P, -theta_zyx, is_real, 'convention', 'zyx', 'rot_intrinsic', rot_intrinsic);

%Block-matrix rotation
D = complex(zeros(size(C)));
for l = 0:P
    idx = ((l^2) + 1):((l+1)^2);
    D(idx, :) = R_list(l+1).R * C(idx, :);
end