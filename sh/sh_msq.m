function D = sh_msq(C, is_real)
%Spherical harmonic expansion magnitude squared expansion
%Y(theta, phi) * D = (Y(theta, phi)*C) * conj(Y(theta, phi)*C) 
%                  = |Y(theta, phi)*C|^2

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:         [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%is_real:   Logical, if true, evaluate real SH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%D:         [(2 * P + 1)^2 x M]  SH coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_msq', '-o', 'sh/sh_msq_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Take the square magnitude of a random function

% rng(541);
% max_odr = 3;
% M = 2;
% C = sh_rand(max_odr, 1);
% D = sh_msq(C);
% sh_plt(C, 'mercator', false, 'dB_lim', [-32, 20]);
% sh_plt(D, 'mercator', false, 'dB_lim', [-32, 20]);

% %Real version
% C = sh_rand(max_odr, 1, true);
% D = sh_msq(C, true);
% sh_plt(C, 'mercator', true, 'dB_lim', [-32, 20]);
% sh_plt(D, 'mercator', true, 'dB_lim', [-32, 20]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare Matlab versus Codegen

% rng(541);
% max_odr = 3;
% M = 2;
% C = sh_rand(max_odr, 1);
% D_mat = sh_msq(C);
% D_mex = sh_msq_mex(C);
% err = norm(D_mat - D_mex)

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
    is_real (1,1) logical = false;
end

if is_real
    C_cpx = sh_re2cpx(C);
    D_cpx = sh_mul(sh_conj(C_cpx), C_cpx);
    D = sh_cpx2re(D_cpx);
    D = real(D);
else
    D = sh_mul(sh_conj(C), C);
end
