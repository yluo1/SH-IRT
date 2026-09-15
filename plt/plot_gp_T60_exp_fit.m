function plot_gp_T60_exp_fit
%Plot examples of exp filter fit to sampling T60 from GP posterior evaluation grid

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup GP prior mean and covariance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options_mu  = gp_mu_opts('mu_func', 'power', 'mu_alpha', 1, 'mu_beta', 0.25, 'mu_fc', 8000);
options_cov = gp_cov_opts('cov_sigma', sqrt(2)/2, 'cov_gamma', 2/3, 'cov_ell', 343 * 1);

Fs = 16000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Azimuth plane spherical coordinate evaluation grid, sampled once, GP posterior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup evaluation grid
N_B = 64;
N_E = 7;
f_lo    = 80;
omega   = 2 * pi * logspace(log10(f_lo), log10(Fs/2), N_B)';
phi     = deg2rad(linspace(0, 180, N_E))';
theta   = pi/2 * ones(N_E, 1);

%Observations
N_S = 3;
omega_obs = logspace(log10(f_lo), log10(Fs/2), N_S)' * 2 * pi;

obs = gp_obs_opts(  'omega', omega_obs, ...
                    'theta', deg2rad(90) * ones(N_S, 1), ...
                    'phi', deg2rad(0) * ones(N_S, 1), ...
                    'T60', mu_lpf(omega_obs, 1, 1, 4000), ...
                    'log_noise_std', 0.005 * ones(N_S, 1) ... %0.5 percent
                    );

options_disp_gp = gp_disp_opts('disp_xlim', [f_lo inf], 'disp_legend_samples', false, ...
    'disp_legend_mean', false, 'disp_legend_var', false, 'disp_legend_loc', 'southwest', ...
    'disp_legend_num_cols', 1, 'disp_ylim', [0, 1.1], ...
    'disp_sample_eval', true, 'disp_mean', false, 'disp_var', true, ...
    'disp_legend_compact', true, 'disp_position', [100, 100, 560 * 0.825, 480 * 0.666], ...
    'disp_legend_transparency', 0.75, 'disp_colororder', 'gem12', ...
    'disp_var_transparency', 0.025, ...
    'disp_sample_stride', 1);

% Sample 1 function
rng(1215 + 10); 
num_evals = 1;

[log_T60, h_fig_posterior_grid] = gp_t60_sample(omega, theta, phi, num_evals, obs, ...
        'options_mu', options_mu, 'options_cov', options_cov, ...
        'enable_disp', true, 'options_disp', options_disp_gp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uniform frequency interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_uni = 32;
omega_uni   = 2 * pi * linspace(0, Fs/2, N_uni)';
T60_uni_sec = exp(interp1(omega, log_T60, omega_uni, 'pchip'));

fontsize = 16;
h_fig_posterior_grid_interp1 = figure; 
semilogx(omega_uni / (2 * pi), T60_uni_sec, 'linewidth', 1.5); 
grid on; axis tight; 
xlabel('Frequency (Hz)', 'fontsize', fontsize); ylabel('T60 (Seconds)', 'fontsize', fontsize); 
title('T60 Interpolation over Uniform Frequencies', 'fontsize', fontsize + 1);
set(gca, 'fontsize', fontsize - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit exponentiating filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_taps = N_uni;
tol0 = 1e-6;

g_RT60  = cell(1, N_E);
g_MP    = cell(1, N_E);
h_g_fig = cell(1, N_E);

idx_fit = N_E;
[g_RT60{idx_fit}, g_MP{idx_fit}, h_g_fig{idx_fit}] = ft_exp_design(N_taps, 'RT60',  'enable_disp', true, 'Fs', Fs, 'RT60_sec', T60_uni_sec(:, idx_fit), 'tol0', tol0);
h_g_fig{idx_fit}.Position = [100, 100, 600, 480];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate noise and filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(521345);
T = 0.5;
h = randn(1, ceil(T * Fs));
[f_1, h_f_exp_fig] = ft_exp_conv_opt(h, g_RT60{N_E}, ...
    'enable_disp', true, 'Fs', Fs, 'N_FFT', 512);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_dir = 'figs/figs_t60';
if ~isfolder(out_dir)
    mkdir(out_dir);
end

exportgraphics(h_fig_posterior_grid, fullfile(out_dir, 'sample_GP_post_exp_fit_field.png'));
exportgraphics(h_fig_posterior_grid_interp1, fullfile(out_dir, 'sample_GP_post_exp_fit_interp.png'));
exportgraphics(h_g_fig{idx_fit}, fullfile(out_dir, 'sample_GP_post_exp_fit_filter.png'));
exportgraphics(h_f_exp_fig, fullfile(out_dir, 'sample_GP_post_exp_fit_exp_conv.png'));