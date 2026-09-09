%function plot_sh_exp_conv_gp_example
%Plot samples of  spherical harmonic exponentiating convolution of expansions C with
%direction-dependent g(theta, phi) sampled from Gaussian process T60(omega, theta, phi)

rng(123);

Fs = 48000;         %Sample rate
T = 1;              %Time (seconds)
M = ceil(Fs * T);   %Number of samples

max_odr = 1;    %Max-order SH expansion of input function
%max_odr_gp = 1; %Max-order SH expansion of directional filtering
%max_odr_gp = 2;
max_odr_gp = 3; %Max-order SH expansion of directional filtering (Default)
%max_odr_gp = 4;

is_real = false;    %C is complex
is_real_dec = true; %Original function is real (pressure field)

%Generate random real function expressed along SH bases
C = sh_rand(max_odr, M, is_real, is_real_dec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify GP mean and coariance priors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options_mu  = gp_mu_opts('mu_func', 'power', 'mu_alpha', 0.5, 'mu_beta', 0); %Constant T60 = 0.5
options_cov = gp_cov_opts('cov_sigma', 0.25, 'cov_gamma', 0.53, 'cov_ell', 343 * 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Sample from prior mean multivariate distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E_mvnrnd_prior = sh_exp_conv_gp(C, [], 'Fs', Fs, 'options_mu', options_mu, 'options_cov', options_cov, ...
%       'max_odr_gp', max_odr_gp, 'sample_method_gp', 'mvnrnd', 'enable_disp', true);
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Sample prior mean function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E_mean_prior = sh_exp_conv_gp(C, [], 'Fs', Fs, 'options_mu', options_mu, 'options_cov', options_cov, ...
%       'max_odr_gp', max_odr_gp, 'sample_method_gp', 'mean', 'enable_disp', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Sample from posterior mean multivariate distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Observations
% N_S = 3;
% omega_obs = logspace(log10(20), log10(Fs/2), N_S)' * 2 * pi;
% 
% obs = gp_obs_opts(  'omega', omega_obs, ...
%                     'theta', deg2rad(90) * ones(N_S, 1), ...
%                     'phi', deg2rad(0) * ones(N_S, 1), ...
%                     'T60', mu_lpf(omega_obs, 2, 1, 4000), ...
%                     'log_noise_std', 0.05 * ones(N_S, 1) ... %5 percent
%                     );

% E_mvnrnd_post = sh_exp_conv_gp(C, obs, 'Fs', Fs, 'options_mu', options_mu, 'options_cov', options_cov, ...
%       'max_odr_gp', max_odr_gp, 'sample_method_gp', 'mvnrnd', 'enable_disp', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Sample prior posterior function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    'disp_legend_mean', false, 'disp_legend_var', false, 'disp_legend_loc', 'southwest', ...
    'disp_legend_num_cols', 1, 'disp_ylim', [0, 1.1], ...
    'disp_sample_eval', true, 'disp_mean', false, 'disp_var', true, ...
    'disp_legend_compact', true, 'disp_position', [100, 100, 560 * 0.825, 480 * 0.666], ...
    'disp_legend_transparency', 0.75, 'disp_colororder', 'gem12', ...
    'disp_var_transparency', 0.025, ...
    'disp_sample_stride', 1);

[E_mean_post, h_fig_gp] = sh_exp_conv_gp(C, obs, 'Fs', Fs, 'options_mu', options_mu, 'options_cov', options_cov, ...
      'max_odr_gp', max_odr_gp, 'sample_method_gp', 'mean', ...
      'N_freq_uni_fit', 24, 'max_taps_g', 24, ...
      'enable_disp', true, 'options_disp_gp', options_disp_gp);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Decode along grid and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_decs = (max_odr + max_odr_gp + 1)^2;
[theta, phi] = sh_fib(N_decs);
[theta_phi_deg] = rad2deg([theta, phi])

f = real(sh_dec(E_mean_post, theta, phi, is_real));

options_plot_RIR = plot_RIR_opts('colormap', hot, 'disp_RIR', false, ...
    'spec_disp_colorbar', true, 'spec_title_interp', 'none', ...
    'spec_disp_yaxis', true, 'fig_size', [400, 600] * (3/4), ...
    'spec_ylim', [50, inf], 'font_size', 18, 'clim', [-80, -20] - 12);


% h_f_dec = cell(1, N_decs);
% for n = 1:N_decs
% %for n = 1:1
%     h_f_dec{n} = plot_RIR(f(n, 1:Fs)', options_plot_RIR);
%     xticks(gca,  [50, 500, 900])
% end
% 
% %Plot for tetrahedron evaluation points
% [theta_tetra, phi_tetra] = sh_plat('tetrahedron');
% f_tetra = real(sh_dec(E_mean_post, theta_tetra, phi_tetra, is_real));
% for n = 1:4
%     plot_RIR(f_tetra(n, 1:Fs)', options_plot_RIR);
% end

%Custom spherical coordinates
theta_cust = deg2rad([90 90 90 90])';
phi_cust = deg2rad([0 60 120 180])';
N_cust = numel(theta_cust);
f_cust = real(sh_dec(E_mean_post, theta_cust, phi_cust, is_real));
h_f_cust_dec = cell(1, N_cust);
for n = 1:N_cust
    h_f_cust_dec{n} = plot_RIR(f_cust(n, 1:ceil(1*Fs))', options_plot_RIR);
    xticks([250 500 750]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_dir = 'figs';
if ~isfolder(out_dir)
    mkdir(out_dir);
end

exportgraphics(h_fig_gp, fullfile(out_dir, 'sample_sh_exp_conv_gp_T60.png'));
for n = 1:N_cust
    exportgraphics(h_f_cust_dec{n} , fullfile(out_dir, ['sample_sh_exp_conv_spec_', num2str(n), '.png']));
end