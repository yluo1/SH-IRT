function tst_sh_fit_msq
%Test sh_fit_msq.m

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:
%tst_sh_fit_msq

%rng(2354); %Original

%N = 250;
N = 81;

max_odr_half = 4;
max_odr = 2 * max_odr_half;

is_real = true;
%is_real = false;

dB_lim = [-60, 0];
disp_phase = false;

is_pdf = true;
%is_pdf = false;

%D_ref = sh_nrm(sh_msq(sh_rand(max_odr_half, 1, is_real), is_real), 'Sum'); %Random function
D_ref = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr_half, pi/2, 2/3, is_real), is_real), 'Sum'); %RBF

[theta, phi] = sh_fib(N);

X = real(sh_dec(D_ref, theta, phi, is_real));

%Solve
[D_fmincon, C_fmincon, err_fmincon] = sh_fit_msq(X, theta, phi, max_odr, is_real, 'fmincon', ...
    'C0', randn((floor(max_odr/2) + 1)^2, 50), 'is_pdf', is_pdf);

[D_sdp, C_sdp, err_sdp] = sh_fit_msq(X, theta, phi, max_odr, is_real, 'sdp', 'is_pdf', is_pdf);

%Plot
sh_plt(D_ref, 'mercator', is_real, 'dB_lim',  dB_lim, 'disp_phase', disp_phase, 'title_name', 'Ref.');

sh_plt(D_fmincon, 'mercator', is_real, 'dB_lim',  dB_lim, 'disp_phase', disp_phase, 'title_name', 'fmincon');
err_fmincon

sh_plt(D_sdp, 'mercator', is_real, 'dB_lim',  dB_lim, 'disp_phase', disp_phase, 'title_name', 'Semi-definite Programming');
err_sdp
