function plot_sh_fit_tr_example
%Plot sh_fit_tr.m modes

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate and sample from random field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(443);

max_odr = 6;
is_real = false;
N_C = (max_odr + 1)^2;
N_pts = (max_odr + 1)^2;

%Reference field
C_ref = sh_rand(max_odr, 1, is_real);

%Sample at spherical coordinates
%[theta, phi] = sh_grd_fib(N_pts);       %Fibonnaci
[theta, phi] = sh_grd_rand(N_pts);  %Random over sphere

X = sh_dec(C_ref, theta, phi, is_real);
%Add noise
X = X + (randn(size(X)) + randn(size(X)) * 1i) * 5e-2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_ls = sh_fit_tr(X, theta, phi, max_odr, is_real, 0);
C_tr_id = sh_fit_tr(X, theta, phi, max_odr, is_real, 0.1, 'identity');
[C_tr_pc, ~, h_tr_pc_list] = sh_fit_tr(X, theta, phi, max_odr, is_real, 0.5, 'picard', 'enable_disp', true);
C_svd_pc = sh_fit_svd(X, theta, phi, max_odr, is_real, 0.5, 'picard');

% Plot reference
dB_lim = [-32, 24];
h_ref = sh_plt(C_ref, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Reference', 'disp_theta_phi', [theta, phi], 'disp_theta_phi_markersize', 12);

h_ls = sh_plt(C_ls, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Non-Regularized Least Squares');

h_tr_id = sh_plt(C_tr_id, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'TR Least Squares L2 Norm');

h_tr_pc = sh_plt(C_tr_pc, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'TR Least Squares Picard');

h_svd_pc = sh_plt(C_svd_pc, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'TSVD Picard');

% Compute errors
err_ls = err_NSHMSQ(C_ref, C_ls)
err_tr_id = err_NSHMSQ(C_ref, C_tr_id)
err_tr_pc = err_NSHMSQ(C_ref, C_tr_pc)
err_svd_pc = err_NSHMSQ(C_ref, C_svd_pc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute errors for varying crossover fractions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_frac = 32;
frac_list = linspace(0, 1, N_frac);

err_tr_pc_list = zeros(1, N_frac);
err_svd_pc_list = zeros(1, N_frac);
for n = 1:N_frac
    C_tr_pc_n = sh_fit_tr(X, theta, phi, max_odr, is_real, frac_list(n), 'picard');
    C_svd_pc_n = sh_fit_svd(X, theta, phi, max_odr, is_real, frac_list(n), 'picard');

    err_tr_pc_list(n) = err_NSHMSQ(C_ref, C_tr_pc_n);
    err_svd_pc_list(n) = err_NSHMSQ(C_ref, C_svd_pc_n);
end

% Plotting

fontsize = 14;
h_pc_frac = figure;
semilogy(frac_list, err_tr_pc_list, 'ro-', frac_list, err_svd_pc_list, 'bd-', 'linewidth', 1.5);
grid on; axis tight;
xlabel('Discrete Picard Condition Crossover Fraction', 'fontsize', fontsize);
ylabel('Normalized Error (NMSQE)', 'fontsize', fontsize);
title('Picard Tikhonov Regularization and Truncated SVD Fitting Error', 'fontsize', fontsize + 1);
set(gca, 'fontsize', fontsize - 1);
h_lg = legend('Tikhonov Regularization', 'Truncated Singular Value Decomposition', 'location', 'best'); set(h_lg, 'fontsize', fontsize - 1);
;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_dir = 'figs/figs_sh/fit';
if ~isfolder(out_dir)
    mkdir(out_dir);
end

% exportgraphics(h_ref{1}, fullfile(out_dir, 'tr_ref.png'));
% exportgraphics(h_tr_pc{1}, fullfile(out_dir, 'tr_pc.png'));
% exportgraphics(h_tr_pc_list{end}, fullfile(out_dir, 'tr_picard.png'));
% exportgraphics(h_pc_frac, fullfile(out_dir, 'tr_tsvd_picard.png'));
% 

