function E = sh_exp_conv(C, g)
%Spherical harmonic exponentiating convolution of expansions C with direction independent g

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:             [(P + 1)^2 x M]   SH coefficients (max order P of M number of functions)
%g:             [1 x N]           Exponentiating FIR filter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%E:             [(P + 1)^2 x M * N]   SH coefficients (max order P of M number of functions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare with exponentiating convolution in pressure domain

% rng(123);
% P = 3;
% M = 1024;
% is_real = true;
% C = sh_rand(P, M, is_real);
% g = ft_two_tap_FIR(0, -1);
% E = sh_exp_conv(C, g);

% %Evaluate at coordinate
% theta = pi / 2; 
% phi = pi /3;
% f_tst = sh_dec(E, theta, phi, is_real);
% f_ref = ft_exp_conv_opt(sh_dec(C, theta, phi, is_real), g);
% err = err_NMSE(f_tst, f_ref)

arguments
    C (:,:) double = complex(0);
    g (1,:) double = 1;
end

N = numel(g);

[N_C, M] = size(C);
P = sqrt(N_C) - 1;
assert(P == floor(P), 'Invalid size C');


E = zeros([N_C, M * N]);

for n = 1:N_C
    E(n, :) = ft_exp_conv_opt(C(n, :), g);
end