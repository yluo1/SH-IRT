function D = sh_cpx2re(C)
%Convert complex spherical harmonic expansion to real spherical harmonic expansion

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:         [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%D:         [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate random complex field, transform to real, evaluate field difference

%rng(123);
%C_cpx = sh_rand(5, 1, false);
%C_re = sh_cpx2re(C_cpx);
%sh_plt(C_cpx, 'mercator', false);
%sh_plt(C_re, 'mercator', true);

%theta = rand(1000,1) * pi;
%phi = rand(1000,1) * 2*pi;
%Y_cpx = sh_dec(C_cpx, theta, phi, false);
%Y_re = sh_dec(C_re, theta, phi, true);
%err = norm(Y_cpx - Y_re)

%err2 = norm(C_cpx - sh_re2cpx(C_re))

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
end

P = size(C, 1);
P = sqrt(P) - 1;
assert(P - floor(P) == 0, 'Invalid size C');

D = complex(zeros(size(C)));

for l = 0:P

    m = (1:l)';
    
    pos = [-1i * ones(l,1); sqrt(2)/2;  (-1).^m] / sqrt(2);    
    neg = [1i * (-1).^flipud(m); sqrt(2)/2; ones(l,1)] / sqrt(2);

    idx = l^2 + (1:(2*l+1));
    D(idx, :) = bsxfun(@times, C(idx, :), pos) + bsxfun(@times, flipud(C(idx, :)), neg);
end
