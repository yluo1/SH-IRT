function plot_gp_optimize_example
% Maximize marginal likelihood w.r.t. hyperparameters and log-T60 samples drawn 
% from another Gaussian process

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup GP prior mean and covariance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options_mu  = gp_mu_opts('mu_func', 'power', 'mu_alpha', 1, 'mu_beta', 0.25, 'mu_fc', 8000);
options_cov = gp_cov_opts('cov_func', 'cov_sqx_chw_ns', 'cov_sigma', 0.5, 'cov_gamma', 0.75, 'cov_ell', 2 * 343);

Fs = 48000; %Sample rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Azimuth plane spherical coordinate evaluation grid, sampled once, GP prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup evaluation grid
N_B     = 16;
N_E     = 24;
omega   = 2 * pi * logspace(log10(20), log10(Fs/2), N_B)';
[theta, phi] = sh_fib(N_E);

options_disp = gp_disp_opts('disp_ylim', [0, 0.6], 'disp_legend_loc', 'southoutside', ...
    'disp_legend_num_cols',  4, 'disp_legend_compact', true, ...
    'disp_position', [100, 100, 800, 800], 'disp_legend_samples', false);

% Sample 1 field
rng(1213); 
num_evals = 1;

[log_T60, h_fig_prior_grid] = gp_t60_sample(omega, theta, phi, num_evals, [], ...
        'options_mu', options_mu, 'options_cov', options_cov, ...
        'enable_disp', true, 'options_disp', options_disp); % [N_B x N_E x num_evals]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify observations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega_grid = repmat(omega, [1, N_E]);
theta_grid = repmat(theta(:)', [N_B, 1]);
phi_grid = repmat(phi(:)', [N_B, 1]);
obs = gp_obs_opts('omega', omega_grid(:), 'theta', theta_grid(:), 'phi', phi_grid(:), 'T60', exp(log_T60(:)), ...
    'log_noise_std', 0.002 * ones(N_B * N_E, 1)); %0.2 percent 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit covariance hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options_fmincon = optimoptions("fmincon", SpecifyObjectiveGradient=true, Display="iter", checkGradients=false, ...
        ScaleProblem=true, ...
        FunctionTolerance=1e-8, ConstraintTolerance=1e-8, OptimalityTolerance=1e-10, StepTolerance=1e-8, ...
        MaxIterations=1000, MaxFunctionEvaluations=10000);

options_cov_0 = gp_cov_opts('cov_func', 'cov_sqx_chw_ns', ...
    'cov_sigma_lim', [0.2, 2], 'cov_gamma_lim', [0.4, 2], 'cov_ell_lim', [100, 1000]);

options_cov_fit = gp_t60_optimize(obs, 'options_mu', options_mu, 'options_cov', options_cov_0, ... 
    'options_fmincon', options_fmincon, 'num_start', 1);

% Compare covariance hyperparameters in structs
disp([newline, 'Prior cov_sigma: ', num2str(options_cov.cov_sigma)]);
disp(['Prior cov_ell: ', num2str(options_cov.cov_ell)]);
disp(['Prior cov_gamma: ', num2str(options_cov.cov_gamma)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample field from posterior distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Draw sample log-T60 from fitted GP
options_disp_fitted = gp_disp_opts('disp_ylim', [0, 0.6], 'disp_legend_loc', 'southoutside', ...
    'disp_legend_num_cols', 4, 'disp_legend_compact', true, ...
    'disp_position', [100, 100, 800, 800], 'disp_legend_samples', false, ...
    'disp_legend_var', false, 'disp_legend_mean', false, ...
    'disp_marker_size', 4);

omega_fitted = 2 * pi * logspace(log10(20), log10(Fs/2), 100)';
[theta_fitted, phi_fitted] = sh_fib(2 * N_E);

[log_T60_fitted, h_fig_post_grid] = gp_t60_sample(omega_fitted, theta_fitted, phi_fitted, num_evals, obs, ...
    'options_mu', options_mu, 'options_cov', options_cov_fit, ...
    'enable_disp', true, 'options_disp', options_disp_fitted);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_dir = 'figs/figs_t60';
if ~isfolder(out_dir)
    mkdir(out_dir);
end

% exportgraphics(h_fig_prior_grid, fullfile(out_dir, 'opt_GP_prior_grid.png'));
% exportgraphics(h_fig_post_grid, fullfile(out_dir,  'opt_GP_post_grid.png'));
