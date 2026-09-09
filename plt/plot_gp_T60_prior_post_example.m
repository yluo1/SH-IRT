function plot_gp_T60_prior_post_example
%Plot examples of sampling T60 from GP prior, posterior, evaluation grids

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setup GP prior mean and covariance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fs = 48000;
options_mu = gp_mu_opts('mu_func', 'power', 'mu_alpha', 1, 'mu_beta', 0.25, 'mu_fc', 8000);
options_cov = gp_cov_opts('cov_sigma', sqrt(2)/2, 'cov_gamma', 2/3, 'cov_ell', 343 * 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Single point spherical coordinate evaluation grid, sampled 4 times, GP prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Setup evaluation grid
N_B = 128;
N_E = 1;
omega = logspace(log10(20), log10(Fs/2), N_B)' * 2 *pi;
[theta, phi] = sh_fib(N_E);
num_evals = 4; %Sample 4 functions

%Sample function
rng(21136 + 8);
[log_T60, h_fig_prior] = gp_t60_sample(omega, theta, phi, num_evals, [], ...
        'options_mu', options_mu, 'options_cov', options_cov, ...
        'enable_disp', true, 'options_disp', gp_disp_opts('disp_ylim', [0, 1.75], 'disp_legend_num_cols', 1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Azimuth plane spherical coordinate evaluation grid, sampled once, GP prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Setup evaluation grid
N_B = 64;
omega = logspace(log10(20), log10(Fs/2), N_B)' * 2 * pi;
phi = deg2rad([0:30:150])';
theta = pi/2 * ones(numel(phi), 1);
num_evals = 1; %Sample 1 function

%Sample function
rng(21136 + 2); 
[log_T60, h_fig_prior_grid] = gp_t60_sample(omega, theta, phi, num_evals, [], ...
        'options_mu', options_mu, 'options_cov', options_cov, ...
        'enable_disp', true, 'options_disp', gp_disp_opts('disp_ylim', [0, 0.75], 'disp_legend_loc', 'northeast', 'disp_legend_num_cols',  2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Single point spherical coordinate evaluation grid, sampled 4 times, GP posterior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Setup evaluation grid
N_B = 128;
N_E = 1;
omega = logspace(log10(20), log10(Fs/2), N_B)' * 2 *pi;
[theta, phi] = sh_fib(N_E);

num_evals = 4; %Sample 4 functions

%Observations
N_S = 3; %3 observations in frequency
omega_obs = logspace(log10(20), log10(Fs/2), N_S)' * 2 * pi;

obs = gp_obs_opts(  'omega', omega_obs, ...
                    'theta', theta(1) * ones(N_S, 1), ...
                    'phi', phi(1) * ones(N_S, 1), ...
                    'T60', mu_pow(omega_obs, 2.5, 0.4), ...
                    'log_noise_std', 0.05 * ones(N_S, 1) ... %5 percent
                    );

%Sample function
rng(21136 + 9);
[log_T60, h_fig_post] = gp_t60_sample(omega, theta, phi, num_evals, obs, ...
            'options_mu', options_mu, 'options_cov', options_cov, ...
            'enable_disp', true, 'options_disp', gp_disp_opts('disp_ylim', [0, 1], 'disp_legend_loc', 'northeast', 'disp_legend_num_cols', 1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Multi-point (2) spherical coordinate evaluation grid, sampled 4 times, GP posterior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Setup evaluation grid
N_B = 128;
N_E = 2;
omega = logspace(log10(20), log10(Fs/2), N_B)' * 2 * pi;
[theta, phi] = sh_fib(N_E);

num_evals = 4; %Sample 4 functions

%Observations
N_S = 3; %3 observations in frequency
omega_obs = logspace(log10(20), log10(Fs/2), N_S)' * 2 * pi;

%Duplicate for two different spherical coordinates
obs = gp_obs_opts(  'omega', repmat(omega_obs, [2, 1]), ...
                    'theta', [theta(1) * ones(N_S, 1); theta(2) * ones(N_S, 1)], ...
                    'phi', [phi(1) * ones(N_S, 1); phi(2) * ones(N_S, 1)], ...
                    'T60', [mu_pow(omega_obs, 2.5, 0.4); mu_pow(omega_obs, 1.5, 0.8) ], ... %Vary the function
                    'log_noise_std', repmat(0.005 * ones(N_S, 1), [2, 1]) ... %0.5 percent
                    );

%Sample function
rng(21136 + 9);
[log_T60, h_fig_post_multi] = gp_t60_sample(omega, theta, phi, num_evals, obs, ...
            'options_mu', options_mu, 'options_cov', options_cov, ...
            'enable_disp', true, 'options_disp', gp_disp_opts('disp_ylim', [0, 1], 'disp_legend_loc', 'southoutside', 'disp_legend_num_cols', 2, 'disp_position', [100, 100, 800, 800]));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_dir = 'figs';
if ~isfolder(out_dir)
    mkdir(out_dir);
end

% exportgraphics(h_fig_prior, fullfile(out_dir, 'sample_GP_prior.png'));
% exportgraphics(h_fig_prior_grid, fullfile(out_dir, 'sample_GP_prior_grid.png'));
% exportgraphics(h_fig_post, fullfile(out_dir, 'sample_GP_post.png'));