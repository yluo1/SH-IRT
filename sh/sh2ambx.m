function D = sh2ambx(C, is_real)
%Convert N3D spherical harmonic expansion to AmbiX format

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:             [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%is_real:       Logical, if true, evaluate real SH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%D:             [(P + 1)^2 x M] AmbiX coefficients (max order P of M number of functions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Check forward and inverse transformations

%is_real = false;
%C_ref = sh_enc_proj(3, pi/2, 0, is_real, true);
%C_amb = sh2ambx(C_ref, is_real);
%C_tst = ambx2sh(C_amb);
%err = norm(C_ref - C_tst)

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
    is_real (1,1) logical = false;
end

C = sh_nrm(C, 'SN3D'); %Schmidt semi-normalization

if is_real
    D = real(C);
else
    D = sh_cpx2re(C);
end