function plot_gp_T60_prior_post_example
%Plot examples of sampling T60 from GP prior, posterior, evaluation grids

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup GP prior mean and covariance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options_mu  = gp_mu_opts('mu_func', 'power', 'mu_alpha', 1, 'mu_beta', 0.25, 'mu_fc', 8000);
options_cov = gp_cov_opts('cov_sigma', sqrt(2)/2, 'cov_gamma', 2/3, 'cov_ell', 343 * 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single point spherical coordinate evaluation grid, sampled 4 times, GP prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs = 48000; %Sample rate

% Setup evaluation grid
N_B = 128; % Number of angular frequencies
N_E = 1;   % Number of spherical coordinates

omega        = 2 * pi * logspace(log10(20), log10(Fs/2), N_B)';
[theta, phi] = sh_fib(N_E);

% Sample 4 functions
num_evals = 4;
rng(21136 + 12); %Seed RNG
[log_T60, h_fig_prior] = gp_t60_sample(omega, theta, phi, num_evals, [], ...
        'options_mu', options_mu, 'options_cov', options_cov, ...
        'enable_disp', true, 'options_disp', gp_disp_opts('disp_ylim', [0, 1.75], 'disp_legend_num_cols', 1));

% for n = 1:20
%     rng(21136 + n);
%     n
%     [log_T60, h_fig_prior] = gp_t60_sample(omega, theta, phi, num_evals, [], ...
%         'options_mu', options_mu, 'options_cov', options_cov, ...
%         'enable_disp', true, 'options_disp', gp_disp_opts('disp_ylim', [0, 1.75], 'disp_legend_num_cols', 1));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Azimuth plane spherical coordinate evaluation grid, sampled once, GP prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup evaluation grid
N_B = 64;
omega   = 2 * pi * logspace(log10(20), log10(Fs/2), N_B)';
phi     = deg2rad([0:30:150])';
theta   = pi/2 * ones(numel(phi), 1);

% Sample 1 function
rng(11136 + 13); 
num_evals = 1;
[log_T60, h_fig_prior_grid] = gp_t60_sample(omega, theta, phi, num_evals, [], ...
        'options_mu', options_mu, 'options_cov', options_cov, ...
        'enable_disp', true, 'options_disp', gp_disp_opts('disp_ylim', [0, 0.75], 'disp_legend_loc', 'northeast', 'disp_legend_num_cols',  2));

% for n = 1:20
%     rng(11136 + n); 
%     [log_T60, h_fig_prior_grid] = gp_t60_sample(omega, theta, phi, num_evals, [], ...
%             'options_mu', options_mu, 'options_cov', options_cov, ...
%             'enable_disp', true, 'options_disp', gp_disp_opts('disp_ylim', [0, 0.75], 'disp_legend_loc', 'northeast', 'disp_legend_num_cols',  2));
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single point spherical coordinate evaluation grid, sampled 4 times, GP posterior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup evaluation grid
N_B = 128;
N_E = 1;
omega        =  2 * pi * logspace(log10(20), log10(Fs/2), N_B)';
[theta, phi] = sh_fib(N_E);

% Observations
N_S = 3; %3 observations in frequency
omega_obs = logspace(log10(20), log10(Fs/2), N_S)' * 2 * pi;

obs = gp_obs_opts(  'omega', omega_obs, ...
                    'theta', theta(1) * ones(N_S, 1), ...
                    'phi', phi(1) * ones(N_S, 1), ...
                    'T60', mu_pow(omega_obs, 2.5, 0.4), ...
                    'log_noise_std', 0.05 * ones(N_S, 1) ... %5 percent
                    );

% Sample function
rng(21136 + 18);
num_evals = 4; %Sample 4 functions
[log_T60, h_fig_post] = gp_t60_sample(omega, theta, phi, num_evals, obs, ...
            'options_mu', options_mu, 'options_cov', options_cov, ...
            'enable_disp', true, 'options_disp', gp_disp_opts('disp_ylim', [0, 1], 'disp_legend_loc', 'northeast', 'disp_legend_num_cols', 1));

% for n = 1:20
%     rng(21136 + n);
%     [log_T60, h_fig_post] = gp_t60_sample(omega, theta, phi, num_evals, obs, ...
%             'options_mu', options_mu, 'options_cov', options_cov, ...
%             'enable_disp', true, 'options_disp', gp_disp_opts('disp_ylim', [0, 1], 'disp_legend_loc', 'northeast', 'disp_legend_num_cols', 1));
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-point (2) spherical coordinate evaluation grid, sampled 4 times, GP posterior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup evaluation grid
N_B = 128;
N_E = 2;
omega = 2 * pi * logspace(log10(20), log10(Fs/2), N_B)';
[theta, phi] = sh_fib(N_E);

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
num_evals = 4; %Sample 4 functions
[log_T60, h_fig_post_multi] = gp_t60_sample(omega, theta, phi, num_evals, obs, ...
            'options_mu', options_mu, 'options_cov', options_cov, ...
            'enable_disp', true, 'options_disp', gp_disp_opts('disp_ylim', [0, 1], 'disp_legend_loc', 'southoutside', 'disp_legend_num_cols', 2, 'disp_position', [100, 100, 800, 800]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Azimuth plane spherical coordinate evaluation grid, sampled once, GP posterior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup evaluation grid
N_B = 64;
omega   = 2 * pi * logspace(log10(20), log10(Fs/2), N_B)';
phi     = deg2rad([0:30:150])';
theta   = pi/2 * ones(numel(phi), 1);

%Observations
N_S = 3;
omega_obs = logspace(log10(20), log10(Fs/2), N_S)' * 2 * pi;

obs = gp_obs_opts(  'omega', omega_obs, ...
                    'theta', deg2rad(90) * ones(N_S, 1), ...
                    'phi', deg2rad(0) * ones(N_S, 1), ...
                    'T60', mu_lpf(omega_obs, 1, 1, 4000), ...
                    'log_noise_std', 0.005 * ones(N_S, 1) ... %0.5 percent
                    );

options_disp_gp = gp_disp_opts('disp_legend_samples', false, ...
    'disp_legend_mean', false, 'disp_legend_var', false, 'disp_legend_loc', 'northeast', ...
    'disp_legend_num_cols', 1, 'disp_ylim', [0, 1.1], ...
    'disp_sample_eval', true, 'disp_mean', false, 'disp_var', true, ...
    'disp_legend_compact', true, 'disp_position', [100, 100, 560 * 0.825, 480 * 0.666], ...
    'disp_legend_transparency', 0.75, 'disp_colororder', 'gem12', ...
    'disp_var_transparency', 0.025, ...
    'disp_sample_stride', 1);

% Sample 1 function
rng(11136 + 13); 
num_evals = 1;
[log_T60, h_fig_posterior_grid] = gp_t60_sample(omega, theta, phi, num_evals, obs, ...
        'options_mu', options_mu, 'options_cov', options_cov, ...
        'enable_disp', true, 'options_disp', options_disp_gp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_dir = 'figs/figs_t60';
if ~isfolder(out_dir)
    mkdir(out_dir);
end

% exportgraphics(h_fig_prior, fullfile(out_dir, 'sample_GP_prior.png'));
% exportgraphics(h_fig_prior_grid, fullfile(out_dir, 'sample_GP_prior_grid.png'));
% exportgraphics(h_fig_post, fullfile(out_dir, 'sample_GP_post.png'));
% exportgraphics(h_fig_posterior_grid, fullfile(out_dir, 'sample_GP_post_grid.png'));