%function plot_sh_ism_example
%Plot examples of spherical harmonic Image-source model RIR augmentation

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate SH-ISM RIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs = 48000;      % Sample rate

max_src_odr = 1; % Max-order source directivity
max_rec_odr = 5; % Max-order receiver directivity

is_real = true;

T = 0.5; % Duration

s = [2 0 1]; % Source coordinates
r = [1 0 1]; % Receiver coordinates
l = [5 6 3]; % Room dimensions

% Room reflection filters with high-frequency dampening
% gamma_pos = [0.8, 0.7, 0.5; 0.2, 0.1, 0.3];
% gamma_neg = [0.9, 0.6, 0.5; 0.1, 0.1, 0.2];

% Room reflection coefficient flat
gamma_pos = db2mag(-[0.8, 1.5, 1.0]);
gamma_neg = db2mag(-[0.5, 2.0, 1.5]);

jitter_coord_bnd = [-1e-1, 1e-1]; % Image coordinate's jitter within +- 10 cm

% Generate SH-ISM
[C, t] = sh_ism(max_src_odr, max_rec_odr, is_real, T, s, r, l, gamma_pos, gamma_neg, ...
    'Fs', Fs, 'jitter_coord_bnd', jitter_coord_bnd);

% Plotting
h_RIR_0 = real(sh_dec(squeeze(C(1, 1, :)).', pi/2, 0, is_real))';
h_RIR_onaxis = real(sh_dec(squeeze(C(1, :, :)), pi/2, 0, is_real))';

h_SH_ism_RIR_0 = plot_RIR(h_RIR_0, plot_RIR_opts('clim', [-120, -60] - 10, 'win_size', 512, 'spec_scale', 'linear'));
h_SH_ism_RIR_onaxis = plot_RIR(h_RIR_onaxis, plot_RIR_opts('clim', [-120, -60], 'win_size', 512, 'spec_scale', 'linear'));
h_SH_ism_SH_hplane = sh_plt(squeeze(C(1, :, :)), 'horizontal', is_real, 'disp_phase', false, 'dB_lim', [-120, 0], 't', t, 'disp_xaxis_ker_size', 1024);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify GP mean and coariance priors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options_mu  = gp_mu_opts('mu_func', 'power', 'mu_alpha', 0.5, 'mu_beta', 0.1);
options_cov = gp_cov_opts('cov_sigma', 0.25, 'cov_gamma', 0.6, 'cov_ell', 343);

options_disp_gp = gp_disp_opts('disp_legend_samples', false, ...
    'disp_legend_mean', false, 'disp_legend_var', false, 'disp_legend_loc', 'southwest', ...
    'disp_legend_num_cols', 1, 'disp_ylim', [0, 0.8], ...
    'disp_sample_eval', true, 'disp_mean', false, 'disp_var', true, ...
    'disp_legend_compact', true, 'disp_position', [100, 100, 560 * 0.825, 480 * 0.666], ...
    'disp_legend_transparency', 0.75, 'disp_colororder', 'gem12', ...
    'disp_var_transparency', 0.025, ...
    'disp_sample_stride', 1);

sample_method_gp = 'mean';

SH_fit_mode = 'svd_ls';
N_freq_uni_fit = 24;
max_taps_g = 24;

% Sample fromr T60 prior mean
obs = [];
[C_sh_ism_T60_prior, h_fig_gp_prior] = sh_exp_conv_gp(squeeze(C(1, :, :)), obs, is_real, ...
    'Fs', Fs, 'options_mu', options_mu, 'options_cov', options_cov, ...
      'max_odr_gp', 3, 'sample_method_gp', sample_method_gp, ...
      'N_freq_uni_fit', N_freq_uni_fit, 'max_taps_g', max_taps_g, ...
      'SH_fit_mode', SH_fit_mode, 'SH_fit_svd_trunc_frac', 0, ...
      'enable_disp', true, 'options_disp_gp', options_disp_gp);

% Sample from T60 posterior mean
N_S = 16;
omega_obs = logspace(log10(20), log10(Fs/2), N_S)' * 2 * pi;
obs = gp_obs_opts(  'omega', omega_obs, ...
                    'theta', deg2rad(90) * ones(N_S, 1), ...
                    'phi', deg2rad(0) * ones(N_S, 1), ...
                    'T60', mu_pow(omega_obs, 0.5, 0), ...
                    'log_noise_std', 0.005 * ones(N_S, 1) ... %0.5 percent
                    );

[C_sh_ism_T60_post, h_fig_gp_post] = sh_exp_conv_gp(squeeze(C(1, :, :)), obs, is_real, 'Fs', Fs, 'options_mu', options_mu, 'options_cov', options_cov, ...
      'max_odr_gp', 3, 'sample_method_gp', sample_method_gp, ...
      'N_freq_uni_fit', N_freq_uni_fit, 'max_taps_g', max_taps_g, ...
      'SH_fit_mode', SH_fit_mode, 'SH_fit_svd_trunc_frac', 0, ...
      'enable_disp', true, 'options_disp_gp', options_disp_gp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decode along grid and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_cust = deg2rad([90 90 90 90])';
phi_cust = deg2rad([0 60 120 180])';
N_cust = numel(theta_cust);

f_cust_orig     = real(sh_dec(squeeze(C(1, :, :)), theta_cust, phi_cust, is_real));
f_cust_prior    = real(sh_dec(C_sh_ism_T60_prior, theta_cust, phi_cust, is_real));
f_cust_post     = real(sh_dec(C_sh_ism_T60_post, theta_cust, phi_cust, is_real));

h_f_cust_dec = cell(3, N_cust);

yticks_list = [250, 1000, 4000, 8000, 16000] / 1000;

for n = 1:N_cust
    options_plot_RIR = plot_RIR_opts('Fs', Fs, 'colormap', hot, 'disp_RIR', false, ...
        'win_size', 64, 'N_FFT', 512, ...
        'spec_disp_colorbar', true, 'spec_title_interp', 'latex', ...
        'spec_disp_yaxis', true, 'fig_size', [600, 300] * (3/4), ...
        'spec_ylim', [100, inf], 'font_size', 18, 'clim', [-160, -50], ...
        'name', ['$(\theta = ',  num2str(rad2deg(theta_cust(n))), '^{\circ}, \phi = ',  num2str(rad2deg(phi_cust(n))), '^{\circ})$'] );

    h_f_cust_dec{1, n} = plot_RIR(f_cust_orig(n, 1:ceil(0.5 * Fs))', options_plot_RIR); yticks(yticks_list);
    h_f_cust_dec{2, n} = plot_RIR(f_cust_prior(n, 1:ceil(0.5 * Fs))', options_plot_RIR); yticks(yticks_list);
    h_f_cust_dec{3, n} = plot_RIR(f_cust_post(n, 1:ceil(0.5 * Fs))', options_plot_RIR); yticks(yticks_list);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_dir = 'figs/figs_t60';
if ~isfolder(out_dir)
    mkdir(out_dir);
end

exportgraphics(h_SH_ism_RIR_0, fullfile(out_dir, 'sh_ism_RIR_0.png'));
exportgraphics(h_SH_ism_RIR_onaxis, fullfile(out_dir, 'sh_ism_RIR_onaxis.png'));
exportgraphics(h_SH_ism_SH_hplane{1}, fullfile(out_dir, 'sh_ism_RIR_hplane.png'));

exportgraphics(h_fig_gp_prior, fullfile(out_dir, 'sh_ism_T60_prior.png'));
exportgraphics(h_fig_gp_post, fullfile(out_dir, 'sh_ism_T60_post.png'));

for n = 1:N_cust
    exportgraphics(h_f_cust_dec{1, n}, fullfile(out_dir, ['sh_ism_orig_', num2str(n), '.png']));
    exportgraphics(h_f_cust_dec{2, n}, fullfile(out_dir, ['sh_ism_T60_prior_', num2str(n), '.png']));
    exportgraphics(h_f_cust_dec{3, n}, fullfile(out_dir, ['sh_ism_T60_post_', num2str(n), '.png']));
end

