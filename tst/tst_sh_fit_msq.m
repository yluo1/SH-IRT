function tst_sh_fit_msq
%Test sh_fit_msq.m

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:
%tst_sh_fit_msq

rng(2354);
N = 100;
max_odr_half = 4;
max_odr = 2 * max_odr_half;
is_real = true;
C_ref = sh_rand(max_odr_half, 1, is_real);
D_ref = sh_msq(C_ref, is_real);
[theta, phi] = sh_fib(N);
X = sh_dec(D_ref, theta, phi, is_real);

K = (max_odr_half + 1)^2;
C0 = randn((floor(max_odr/2) + 1)^2, K);
B0 = randn(K, 2 * K);


[D_MS, C_MS, ~]                     = sh_fit_msq(X, theta, phi, max_odr, is_real, 'MS', 'C0', [C0]);
[D_SOMS, C_SOMS, ~]                 = sh_fit_msq(X, theta, phi, max_odr, is_real, 'SOMS');
[D_MOMS, C_MOMS, ~, output_MOMS]    = sh_fit_msq(X, theta, phi, max_odr, is_real, 'MOMS', 'B', [C0]);
[D_MP, C_MP, ~, output_MP]          = sh_fit_msq(X, theta, phi, max_odr, is_real, 'MP', 'B', [C0], 'B0', B0);

err_MS = err_SHMSQ(D_ref, D_MS)
err_SOMS = err_SHMSQ(D_ref, D_SOMS)
err_MOMS = err_SHMSQ(D_ref, D_MOMS)
err_MP = err_SHMSQ(D_ref, D_MP)

dB_lim = [-80, 20];

% sh_plt(C_ref, 'mercator', is_real, 'title_name', 'Ref. C', 'disp_theta_phi', [theta, phi], 'dB_lim', dB_lim);
% sh_plt(C_MS, 'mercator', is_real, 'title_name', 'Magnitude Square C', 'disp_theta_phi', [theta, phi], 'dB_lim', dB_lim);
% sh_plt(C_SOMS, 'mercator', is_real, 'title_name', 'Sum-of-Magnitude Square C', 'disp_theta_phi', [theta, phi], 'dB_lim', dB_lim);

sh_plt(D_ref, 'mercator', is_real, 'title_name', 'Ref. D', 'disp_theta_phi', [theta, phi], 'dB_lim', dB_lim);
sh_plt(D_MS, 'mercator', is_real, 'title_name', 'Magnitude Square D', 'disp_theta_phi', [theta, phi], 'dB_lim', dB_lim);
sh_plt(D_SOMS, 'mercator', is_real, 'title_name', 'Sum-of-Magnitude Square  D', 'disp_theta_phi', [theta, phi], 'dB_lim', dB_lim);
sh_plt(D_MOMS, 'mercator', is_real, 'title_name', 'Mix-of-Magnitude Square D', 'disp_theta_phi', [theta, phi], 'dB_lim', dB_lim);
sh_plt(D_MP, 'mercator', is_real, 'title_name', 'Mixture Power D', 'disp_theta_phi', [theta, phi], 'dB_lim', dB_lim);
