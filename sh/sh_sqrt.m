function [D, E, F] = sh_sqrt(C, is_real, mode)
%Spherical harmonic expansion square root
%Y(theta, phi) * D \approx sqrt(Y(theta, phi) * C)

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:         [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%is_real:   Logical, if true, evaluate real SH
%mode:      String, method {'sqrt', 'mag_sqrt'}
%                           'sqrt':         Y(theta, phi) * D \approx sqrt(Y(theta, phi) * C)
%                           'mag_sqrt':     Y(theta, phi) * D \approx sqrt(abs(Y(theta, phi) * C))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%D:         [(2 * P + 1)^2 x M]         SH coefficients
%E:         [(P + 1)^2 x M]             SH coefficients fitted to max-order P
%F:         [(floor(P/2) + 1)^2 x M]    SH coefficients fitted to max-order floor(P/2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_sqrt', '-o', 'sh/sh_sqrt_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare magnitude square and square root of squared SH expansions

%is_real = false;
%dB_lim = [-60, 0];
%C = sh_enc_rbf('SqExp', 3, pi/3, pi/3, 0.2, is_real);
%C_sq = sh_msq(C, is_real);
%[D_sq_sqrt, E_sq_sqrt, F_sq_sqrt] = sh_sqrt(C_sq, is_real, 'sqrt');

%err_sqrt = norm([C; zeros(size(D_sq_sqrt, 1) - size(C, 1), 1)] - D_sq_sqrt)

%sh_plt(C, 'mercator', is_real, 'dB_lim', dB_lim)
%sh_plt(D_sq_sqrt, 'mercator', is_real, 'dB_lim', dB_lim)
%sh_plt(E_sq_sqrt, 'mercator', is_real, 'dB_lim', dB_lim)
%sh_plt(F_sq_sqrt, 'mercator', is_real, 'dB_lim', dB_lim)

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
    is_real (1,1) logical = false;
    mode (1,:) char {mustBeMember(mode, {'sqrt', 'mag_sqrt'})} = 'sqrt';
end

[P, M] = size(C);
P = sqrt(P) - 1;
assert(P - floor(P) == 0, 'Invalid size C');

%Fit max-order 2*P
[theta, phi] = sh_fib((2 * P + 1)^2);
if strcmp(mode, 'sqrt')
    f = sqrt(sh_dec(C, theta, phi, is_real));
elseif strcmp(mode, 'mag_sqrt')
    f = sqrt(abs(sh_dec(C, theta, phi, is_real)));
else
    error('Unsupported mode');
end
Y = sh_val(2 * P, theta, phi, is_real);
D = Y \ f;

if nargout > 1
    %Fit max-order P
    E = Y(:, 1:(P + 1)^2) \ f;
end

if nargout > 2
    %Fit max-order floor(P/2)
    F = Y(:, 1:(floor(P / 2) + 1)^2) \ f;
end

