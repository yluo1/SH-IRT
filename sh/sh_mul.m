function [E, A] = sh_mul(C, D)
%Spherical harmonic expansion multiplication via Clebsch-Gordan coefficients

%f(theta, phi) = C.' * Y(theta, phi)
%g(theta, phi) = D.' * Y(theta, phi)

%f(theta, phi) * g(theta, phi) = E.' * Y(theta, phi)

%Express output SH coefficients in terms of transfer matrix A vector product:
%E(:, m) = A(:, :, m) * D(:, m) forall 1 <= m <= M

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:             [(P_C + 1)^2 x M]   SH coefficients (max order P_C of M number of functions)
%D:             [(P_D + 1)^2 x M]   SH coefficients (max order P_D of M number of functions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ouput
%E:             [(P_C + P_D + 1)^2 x M]                     SH coefficients
%A:             [(P_C + P_D + 1)^2 x (P_D + 1)^2 x M]       Pre-computed Clebsch-Gordan \diamond C matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_mul', '-o', 'sh/sh_mul_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare Matlab versus Mex

%rng(2134);
%P_C = 3;
%P_D = 5;
%M = 4;
%C = sh_rand(P_C, M);
%D = sh_rand(P_D, M);

%tic; [E_mat, A_mat] = sh_mul(C, D); toc_mat = toc
%tic; [E_mex, A_mex] = sh_mul_mex(C, D); toc_mex = toc

%err_E = norm(E_mat - E_mex)
%err_A = isequal(A_mat, A_mex)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Spatial filter random SH expansion by radial basis function

% rng(2134);
% P_C = 3;
% P_D = 3;
% theta = pi/2;
% phi = 0;
% ell = 1;
% C = sh_enc_rbf('SqExp', P_C, theta, phi, ell, false, false);
% D = sh_rand(P_D, 1);
% [E, A] = sh_mul(C, D);
% sh_plt(C, 'mercator', false, 'title_name', 'Spatial Filter');
% sh_plt(D, 'mercator', false, 'title_name', 'Random');
% sh_plt(E, 'mercator', false, 'title_name', 'Random * Directional Filter');

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0)
    D (:,:) double {coder.mustBeComplex} = complex(0)
end

[P_C, M_C] = size(C);
[P_D, M_D] = size(D);
P_C = sqrt(P_C) - 1;
P_D = sqrt(P_D) - 1;

assert(M_C == M_D, 'M not matched for C, D');
assert(P_C - floor(P_C) == 0, 'Invalid size C');
assert(P_D - floor(P_D) == 0, 'Invalid size D');

A = complex(zeros( [(P_C + P_D + 1)^2, (P_D + 1)^2, M_C] )); %Pre-compute transfer matrix: Clebsch-Gordan \diamond C matrix

for L = 0:(P_C + P_D)
    for M = -L : L
        idx_A = (L^2) + M + L + 1;

        for l = 0:P_C
             for k = abs(L - l) : min(abs(L + l), P_D)     
                 c = sqrt((2 * l + 1) * (2 * k + 1)/ (4 * pi)) * (-1)^(M) * sqrt(2 * L + 1);
                 cw = c * Wigner3j(l, k, L, 0, 0, 0);

                 for m = max(M - k, -l) : min(M + k, l)
                    cww = cw * Wigner3j(l, k, L, m, M - m, M);

                    if cww ~= 0
                        idx_C   = (l^2) + m + l + 1;
                        idx_A_n = (k^2) + M - m + k + 1;
                        A(idx_A, idx_A_n, :) = A(idx_A, idx_A_n, :)  + cww  *  reshape(C(idx_C, :), [1, 1, M_C]);
                    end
                end
            end
        end

    end
end

E = complex(zeros((P_C + P_D + 1)^2, M_C));
for m = 1:M_C
    E(:, m) = A(:, :, m) * D(:, m);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wigner3j.m by David Terr, Raytheon, 6-17-04
% Compute the Wigner 3j symbol using the Racah formula [1]. 
% Reference: Wigner 3j-Symbol entry of Eric Weinstein's Mathworld: http://mathworld.wolfram.com/Wigner3j-Symbol.html

function wigner = Wigner3j(j1, j2, j3, m1, m2, m3)

t1 = j2 - m1 - j3;
t2 = j1 + m2 - j3;
t3 = j1 + j2 - j3;
t4 = j1 - m1;
t5 = j2 + m2;
tmin = max( 0, max( t1, t2 ) );
tmax = min( t3, min( t4, t5 ) );
wigner = 0;
for t = tmin:tmax
    wigner = wigner + (-1)^t / ( factorial(t) * factorial(t-t1) * factorial(t-t2) ...
        * factorial(t3-t) * factorial(t4-t) * factorial(t5-t) );
end
wigner = wigner * (-1)^(j1-j2-m3) ...
    * sqrt( factorial(j1+j2-j3) * factorial(j1-j2+j3) * factorial(-j1+j2+j3) / factorial(j1+j2+j3+1)...
        * factorial(j1+m1) * factorial(j1-m1) * factorial(j2+m2) * factorial(j2-m2) * factorial(j3+m3) * factorial(j3-m3) );
