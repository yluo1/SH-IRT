function tst_com_sh_dec_rec_conv
%Test commutative properties of spherical harmonic expansion decoding and 
%recursive convolution

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(623);
Fs = 48000;

%T = 1;
%T = 0.1;
T = 0.01;
M = ceil(T * Fs);

max_odr = 5;
N = (max_odr + 1)^2;

is_real = true;
C = sh_rand(max_odr, M, is_real);
g = ft_two_tap_FIR(0, -6, false);

theta = deg2rad(33);
phi = deg2rad(50);

h = sh_dec(C, theta, phi, is_real);

f_pre = ft_rec_conv_opt(h, g); %Pre-decode

D = zeros(N, M * numel(g));
for n = 1:N
    D(n, :) = ft_rec_conv_opt(C(n,:), g);
end
f_pst = sh_dec(D, theta, phi, is_real); %Post-decode


figure;
t = 1:numel(f_pre);
plot(t, f_pre, t, f_pst, t, f_pre - f_pst, 'LineWidth', 1.5);
grid on; axis tight;
legend({'pre', 'pst', 'pre - pst'}, 'Location', 'best');

err = err_SNMSE(f_pre, f_pst) %Does commute