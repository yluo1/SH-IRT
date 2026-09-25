function is_valid_pdf = sh_pdf_check(C_pdf, is_real, tol0)
%Partial check if spherical harmonic expansion is valid probability density
%SH expands into a real field, integrates to unity.

%Note: Omits check on non-negative density, as it requires expensive
%sum-of-magnitude squared factorization of the SH expansion

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:                 [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%is_real:           Logical, if true, evaluate real SH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%is_valid_pdf:      Logical, true if all columns of C_pdf have SH expansions that are non-negative and sum to unity 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Test sample SH RBF expansion

%max_odr = 3;
%theta = pi/2;
%phi = 0;
%ell = 1;
%is_real = false;

%C = sh_enc_rbf('SqExp', max_odr, theta, phi, ell, is_real);
%is_valid_pdf = sh_pdf_check(C)

%C = sh_msq(C, is_real);
%is_valid_pdf = sh_pdf_check(C)

%C = sh_nrm(C, 'Sum');
%is_valid_pdf = sh_pdf_check(C)

arguments
     C_pdf (:,:) double {coder.mustBeComplex} = complex(0);
     is_real (1,1) logical = false;
     tol0 (1,1) double {mustBeNonnegative} = 1e-6;
end

%Check if field is real
if is_real
    C_pdf_cpx = sh_re2cpx(C_pdf);
    is_real_field = all(vecnorm(C_pdf_cpx - sh_conj(C_pdf_cpx)) <= tol0);
else
    is_real_field = all(vecnorm(C_pdf - sh_conj(C_pdf)) <= tol0);    
end
%Note: non-negative requires sum-of-square factorization of expansion

%Check if field integrates to unity
is_unity_integral = all(abs(sh_int(C_pdf, 'Sum') - 1) <= tol0);

is_valid_pdf = is_real_field && is_unity_integral;
