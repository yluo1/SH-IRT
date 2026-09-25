function err = err_NSHMSQ(C, D)
%Spherical harmonic mean magnitude squared error normalized by variance of C over spherical coordinates

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:             [(P_C + 1)^2 x M]   SH coefficients (max order P_C of M number of functions)
%D:             [(P_D + 1)^2 x M]   SH coefficients (max order P_D of M number of functions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%err:           [1 x M]

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
    D (:,:) double {coder.mustBeComplex} = complex(0);    
end

[N_C, M_C] = size(C);
[N_D, M_D] = size(D);

P_C = sqrt(N_C) - 1;
P_D = sqrt(N_D) - 1;

assert(M_C == M_D, 'Size M mismatch in C, D');
assert(P_C - floor(P_C) == 0, 'Invalid size C');
assert(P_D - floor(P_D) == 0, 'Invalid size D');

P = max(P_C, P_D); % Max-order 

E = sh_resize(C, P) - sh_resize(D, P);

err = sh_int(E, 'PowAvg') ./ sh_int(C, 'Var');