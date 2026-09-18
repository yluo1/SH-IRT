function plot_sh_rand_pp_example
%Plot examples of spherical harmonic Poisson-process RIR augmentation

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate probability density functions over spherical coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_odr = 12;
is_real = true;

T	= 0.5;          % Duration
Fs	= 16000;        % Sample rate
M	= ceil(T * Fs); % Number of samples

% Squared exponential of chordal distance radial basis function
C_pdf = sh_nrm(sh_enc_rbf('SqExp', max_odr, pi/2, 0, 0.25, is_real), 'Sum');

C_pdf = C_pdf * ones(1, M);
D_lin = sh_pdf_scatter(C_pdf, 'lin', 0.999, is_real);
D_exp = sh_pdf_scatter(C_pdf, 'exp', 0.999, is_real, 'exp_k', 0.1);

dB_lim = [-40, 20];
t = (0:(M-1)) / Fs;
h_pdf = sh_plt(C_pdf, 'horizontal', is_real, 'disp_phase', false, 'title_name_override', 'Squared Exponential Kernel Density', 'dB_lim', dB_lim, 't', t);
h_pdf_lin = sh_plt(D_lin, 'horizontal', is_real, 'disp_phase', false, 'title_name_override', 'Linear Scattering Density', 'dB_lim', dB_lim, 't', t);
h_pdf_exp = sh_plt(D_exp, 'horizontal', is_real, 'disp_phase', false, 'title_name_override', 'Exponential k = 0.5 Scattering Density', 'dB_lim', dB_lim, 't', t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate absolute echo density profile and sample SRIRs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aed = logspace(log10(10/Fs), log10(2), M) * Fs;
rseed = 12543 + 7;

rng(rseed);
[C_pp, t, h_pp] = sh_rand_pp(max_odr, aed, is_real, 'Fs', Fs, ...
    'direct_unity_first_pulse', true, 'C_pdf', C_pdf, 'enable_disp', true);

rng(rseed);
[C_pp_lin, t, h_pp_lin] = sh_rand_pp(max_odr, aed, is_real, 'Fs', Fs, ...
    'direct_unity_first_pulse', true, 'C_pdf', D_lin, 'enable_disp', true);

rng(rseed);
[C_pp_exp, t, h_pp_exp] = sh_rand_pp(max_odr, aed, is_real, 'Fs', Fs, ...
    'direct_unity_first_pulse', true, 'C_pdf', D_exp, 'enable_disp', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify GP mean and coariance priors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options_mu  = gp_mu_opts('mu_func', 'power', 'mu_alpha', 0.5, 'mu_beta', 0.1);
options_cov = gp_cov_opts('cov_sigma', 0.25, 'cov_gamma', 0.75, 'cov_ell', 343);

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
[C_pp_exp_T60_prior, h_fig_gp_prior] = sh_exp_conv_gp(C_pp_exp, obs, is_real, ...
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

[C_pp_exp_T60_post, h_fig_gp_post] = sh_exp_conv_gp(C_pp_exp, obs, is_real, 'Fs', Fs, 'options_mu', options_mu, 'options_cov', options_cov, ...
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

f_cust_orig     = real(sh_dec(C_pp_exp, theta_cust, phi_cust, is_real));
f_cust_prior    = real(sh_dec(C_pp_exp_T60_prior, theta_cust, phi_cust, is_real));
f_cust_post     = real(sh_dec(C_pp_exp_T60_post, theta_cust, phi_cust, is_real));

h_f_cust_dec = cell(3, N_cust);

yticks_list = [250, 1000, 4000] / 1000;

for n = 1:N_cust
    options_plot_RIR = plot_RIR_opts('Fs', Fs, 'colormap', hot, 'disp_RIR', false, ...
        'win_size', 64, 'N_FFT', 512, ...
        'spec_disp_colorbar', true, 'spec_title_interp', 'latex', ...
        'spec_disp_yaxis', true, 'fig_size', [600, 300] * (3/4), ...
        'spec_ylim', [50, inf], 'font_size', 18, 'clim', [-160, -50], ...
        'name', ['$(\theta = ',  num2str(rad2deg(theta_cust(n))), '^{\circ}, \phi = ',  num2str(rad2deg(phi_cust(n))), '^{\circ})$'] );

    h_f_cust_dec{1, n} = plot_RIR(f_cust_orig(n, 1:ceil(0.5 * Fs))', options_plot_RIR);  yticks(yticks_list);
    h_f_cust_dec{2, n} = plot_RIR(f_cust_prior(n, 1:ceil(0.5 * Fs))', options_plot_RIR);  yticks(yticks_list);
    h_f_cust_dec{3, n} = plot_RIR(f_cust_post(n, 1:ceil(0.5 * Fs))', options_plot_RIR);  yticks(yticks_list);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_dir = 'figs/figs_t60';
if ~isfolder(out_dir)
    mkdir(out_dir);
end

% exportgraphics(h_pdf{1}, fullfile(out_dir, 'rand_pp_pdf.png'));
% exportgraphics(h_pdf_lin{1}, fullfile(out_dir, 'rand_pp_pdf_lin.png'));
% exportgraphics(h_pdf_exp{1}, fullfile(out_dir, 'rand_pp_pdf_exp.png'));
% 
% exportgraphics(h_pp{1}, fullfile(out_dir, 'SRIR_pp_pdf.png'));
% exportgraphics(h_pp{2}, fullfile(out_dir, 'SRIR_hplane_pp_pdf.png'));
% 
% exportgraphics(h_pp_lin{1}, fullfile(out_dir, 'SRIR_pp_pdf_lin.png'));
% exportgraphics(h_pp_lin{2}, fullfile(out_dir, 'SRIR_hplane_pp_pdf_lin.png'));
% 
% exportgraphics(h_pp_exp{1}, fullfile(out_dir, 'SRIR_pp_pdf_exp.png'));
% exportgraphics(h_pp_exp{2}, fullfile(out_dir, 'SRIR_hplane_pp_pdf_exp.png'));

% exportgraphics(h_fig_gp_prior, fullfile(out_dir, 'SRIR_pp_gp_exp_T60_prior.png'));
% exportgraphics(h_fig_gp_post, fullfile(out_dir, 'SRIR_pp_gp_exp_T60_post.png'));
% 
for n = 1:N_cust
    exportgraphics(h_f_cust_dec{1, n}, fullfile(out_dir, ['SRIR_pp_gp_exp_orig_', num2str(n), '.png']));
    exportgraphics(h_f_cust_dec{2, n}, fullfile(out_dir, ['SRIR_pp_gp_exp_T60_prior_', num2str(n), '.png']));
    exportgraphics(h_f_cust_dec{3, n}, fullfile(out_dir, ['SRIR_pp_gp_exp_T60_post_', num2str(n), '.png']));
end
