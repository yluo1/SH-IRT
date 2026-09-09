function D = sh_conj(C)
%Conjugate spherical harmonic expansion of function f(theta, phi):

%f(theta, phi)          = C.' * Y(theta, phi)
%conj(F(theta, phi))    = D.' * Y(theta, phi)

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:         [(P+1)^2 x M] SH coefficients (max order P of M number of functions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%D:         [(P+1)^2 x M] SH coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation

%codegen('sh_conj', '-o', 'sh/sh_conj_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Evaluate function and its conjugate after decoding

% rng(534);
% P = 4;
% M = 5;
% N = 10;
% theta = rand(N, 1) * pi;
% phi = rand(N, 1) * 2 * pi;
% C = sh_rand(P, M);
% D = sh_conj(C);
% f = sh_dec(C, theta, phi);
% f_conj = sh_dec(D, theta, phi);
% err = norm(conj(f) - f_conj)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare Matlab versus Mex

% rng(534);
% P = 4;
% M = 5;
% C = sh_rand(P, M);
% D_mat = sh_conj(C);
% D_mex = sh_conj(C);
% err = norm(D_mat - D_mex)

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
end

[P, M] = size(C);
P = sqrt(P) - 1;
assert(P - floor(P) == 0, 'Invalid size C');

D = complex(zeros((P+1)^2, M));
for l = 0:P
    idx_D = (l^2 + 1) : (l + 1)^2;
    idx_C = fliplr(idx_D);

    D(idx_D, :) = bsxfun(@times, conj(C(idx_C, :)), ((-1).^(-l:l))' );
end