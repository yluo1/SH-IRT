function D = sh_resize(C, max_odr)
%Resize max-odr of spherical harmonic expansions

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:             [(P+1)^2 x M] SH coefficients (max order P of M number of functions)
%max_odr:       Max order expansion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%D:             [(max_odr+1)^2 x M] SH coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Truncate and zero-pad expansions

% rng(534);
% P = 4;
% M = 1;
% is_real = false;
% C = sh_rand(P, M, is_real);
% D_2 = sh_resize(C, 2);
% D_6 = sh_resize(C, 6);

% sh_plt(C, 'mercator', is_real);
% sh_plt(D_2, 'mercator', is_real);
% sh_plt(D_6, 'mercator', is_real);


arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 3;
end

[P, M] = size(C);
P = sqrt(P) - 1;
assert(P - floor(P) == 0, 'Invalid size C');

D = [ C(1:(min(P, max_odr) + 1)^2, :); zeros([(max_odr+1)^2 - (P+1)^2, M]) ];


