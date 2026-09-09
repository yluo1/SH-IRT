function E = sh_conv(C, D)
%Spherical harmonic convolution of expansions C with D

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:             [(P + 1)^2 x M]   SH coefficients (max order P of M number of functions)
%D:             [(P + 1)^2 x M]   SH coefficients (max order P of M number of functions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ouput
%E:             [(P + 1)^2 x M]                     SH coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_conv', '-o', 'sh/sh_conv_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Blurring random function with RBF in SH domain

% rng(5451)
% max_odr = 16;
% theta = 0;
% phi = 0;
% ell = 0.15;
% is_real = false;
% 
% C_rand = sh_rand(max_odr, 1, is_real);
% C_RBF = sh_enc_rbf('SqExpNorm', max_odr, theta, phi, ell, is_real, false);
% C_conv = sh_conv(C_rand, C_RBF);
% dB_lim = [-24, 24];
% sh_plt(C_RBF, 'mercator', false, 'dB_lim', dB_lim);
% sh_plt(C_rand, 'mercator', false, 'dB_lim', dB_lim);
% sh_plt(C_conv, 'mercator', false, 'dB_lim', dB_lim);
% var_rand = sh_int(C_rand, 'Var')
% var_conv = sh_int(C_conv, 'Var') %Reduced variance

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
    D (:,:) double {coder.mustBeComplex} = complex(0);
end

[P_C, M_C] = size(C);
[P_D, M_D] = size(D);
P_C = sqrt(P_C) - 1;
P_D = sqrt(P_D) - 1;

assert(P_C == P_D && M_C == M_D, 'Size C, D mismatch');
assert(P_C - floor(P_C) == 0, 'Invalid size C');

P = P_C;
M = M_C;
E = complex((P + 1)^2, M);
for l = 0:P_C
    idx_C = (l^2 + 1) : (l+1)^2;
    idx_D = l^2 + l + 1;
    E(idx_C, :) = sqrt(4 * pi / (2 * l + 1) ) * bsxfun(@times, C(idx_C, :), D(idx_D, :));
end

%Note: Expressible as diagonal matrix-vector product
