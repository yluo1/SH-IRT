function tst_dir_dep_exp_conv
%Test direction-dependent exponentiating convolution SH encoding and decoding

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(2343);
Fs = 48000;

%T = 1;
T = 0.5;
%T = 0.1;
%T = 0.01;

M = ceil(T * Fs); 
max_odr = 3;
is_real = true;

%Random pressure field SH expansion
C = sh_rand(max_odr, M, is_real);

%Direction-dependent two-tap exponentiating filter
theta_RT60_DC = deg2rad(90);
phi_RT60_DC = deg2rad(77);
ell_RT60_DC = 1;
max_RT60_DC = 0.5;

theta_RT60_NQ = deg2rad(90);
phi_RT60_NQ = deg2rad(77);
ell_RT60_NQ = 1;
max_RT60_NQ = 0.25;

N_g = 2; %Two tap exponentiating FIR

%Max-order depends on SH expansion order of the direction-dependent RT60 functions under exponentiation
%Minimally should be max_odr(C) + max_odr(D)

%max_odr_grid = 1 * max_odr;  
%max_odr_grid = 2 * max_odr;
max_odr_grid = 3 * max_odr;

N_grid = (max_odr_grid + 1)^2; %Minimum number of points on grid

[theta_grid, phi_grid] = sh_fib(N_grid); %Evaluation grid

RT60_DC_grid = max_RT60_DC * srbf_val('SqExp', theta_RT60_DC, phi_RT60_DC, ell_RT60_DC, theta_grid, phi_grid); %[N_grid x 1]
RT60_NQ_grid = max_RT60_NQ * srbf_val('SqExp', theta_RT60_NQ, phi_RT60_NQ, ell_RT60_NQ, theta_grid, phi_grid); %[N_grid x 1]
h_grid       = sh_dec(C, theta_grid, phi_grid, is_real); %[N_grid x M]
hg_grid      = zeros(N_grid, M * N_g);

%Compute time-varying convolution at each grid point
for n = 1:N_grid
    g_n = ft_exp_design(2, 'RT60', 'RT60_sec', [RT60_DC_grid(n), RT60_NQ_grid(n)], 'fit_method', 'two_tap');
    hg_grid(n, :) = ft_exp_conv_opt(h_grid(n, :), g_n);
end

%CG = sh_fit_svd(hg_grid, theta_grid, phi_grid, max_odr_grid, is_real); 
CG = sh_val(max_odr_grid, theta_grid, phi_grid, is_real) \ hg_grid;     %Projection onto SH


%Compute solutions at evaluation point
% theta_eval = deg2rad(70);
% phi_eval = deg2rad(60);

theta_eval = deg2rad(50);
phi_eval = deg2rad(80);

RT60_DC_eval = max_RT60_DC * srbf_val('SqExp', theta_RT60_DC, phi_RT60_DC, ell_RT60_DC, theta_eval, phi_eval); 
RT60_NQ_eval = max_RT60_NQ * srbf_val('SqExp', theta_RT60_NQ, phi_RT60_NQ, ell_RT60_NQ, theta_eval, phi_eval);
g_eval          = ft_exp_design(2, 'RT60', 'RT60_sec', [RT60_DC_eval, RT60_NQ_eval], 'fit_method', 'two_tap');
h_eval          = sh_dec(C, theta_eval, phi_eval, is_real);


f_ref = ft_exp_conv_opt(h_eval, g_eval);                %Reference
f_tst = sh_dec(CG, theta_eval, phi_eval, is_real);      %Observed

err = err_NMSE(f_tst, f_ref)

% plot_RIR(f_ref);
% plot_RIR(f_tst);

%Error plotting
figure; t = 1:(M * N_g); plot(t, f_ref, t, f_tst, t, f_ref - f_tst, 'k', 'linewidth', 1.5);
grid on; axis tight;
legend('ref', 'tst', 'ref - tst');


