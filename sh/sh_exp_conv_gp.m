function [E, h_fig_gp] = sh_exp_conv_gp(C, obs, is_real, options)
%Spherical harmonic exponentiating convolution of expansions C with
%direction-dependent exponentiating filter g(theta, phi) sampled from
%Gaussian process T60(omega, theta, phi)

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:             [(P + 1)^2 x M]   SH coefficients (max order P of M number of functions)

%obs:                   Struct, observed T60(omega, theta, phi), see gp_obs_opts.m
%                       [] for none

%is_real:               Logical, if true, evaluate real SH

%options:                  Struct

%options.Fs:               Sampling rate

%options.options_mu_gp:    Struct, prior mean options, see gp_mu_opts.m
%options.options_cov_gp:   Struct, prior covariance options, see gp_cov_opts.m

%options.max_odr_gp:        Max SH expansion order for functions drawn from GP
%options.sample_method_gp:  String, GP sampling method {'mvnrnd', 'mean'}
%options.N_freq_val_gp:     Number of uniform log-frequency points from DC to Nyquist to sample from GP

%options.N_freq_uni_fit:    Number of uniform frequency points to interpolate sampled T60
%options.max_taps_g:        Number of filter taps to fit exponentiating filter per T60

%options.SH_fit_mode:               String, fitting method {'svd_inv', 'svd_ls'}
%                                       'svd_inv':      Truncated inverse
%                                       'svd_ls'        Truncated leasts-squares
%options.SH_fit_svd_trunc_frac:     SH fit truncates fraction of largest singular values

%options.enable_disp:       Logical, if true, plot GP prior or posterior

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%E:  [(P + options.max_odr_gp + 1)^2 x (M * max_taps_g)]   SH coefficients for 'svd_inv' options.SH_fit_mode
%    [(P + 1)^2 x (M * max_taps_g)]                        SH coefficients for 'svd_ls'  options.SH_fit_mode

%h_fig_gp:  Handle to GP figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:
% rng(123);
% Fs = 48000;
% max_odr = 1;
% max_odr_gp = 1;
% M = Fs * 1;
% is_real = false;
% is_real_dec = true;
% C = sh_rand(max_odr, M, is_real, is_real_dec);
% obs = [];
% options_mu  = gp_mu_opts('mu_func', 'LPF', 'mu_beta', 1);
% options_cov = gp_cov_opts('cov_sigma', 0.1, 'cov_gamma', 2/3);
% %sample_method_gp = 'mvnrnd';
% sample_method_gp = 'mean';
% E = sh_exp_conv_gp(C, obs, is_real, 'Fs', Fs, 'options_mu', options_mu, 'options_cov', options_cov, ...
%       'max_odr_gp', max_odr_gp, 'sample_method_gp', sample_method_gp, 'enable_disp', true);

% %Decode
% [theta, phi] = sh_plat('tetrahedron');
% f = real(sh_dec(E, theta, phi, is_real));
% plot_RIR(f(1, 1:Fs)');
% plot_RIR(f(2, 1:Fs)');

arguments
    C (:,:) double = 0;

    obs = [];

    is_real (1,1) logical = false;

    options.Fs (1,1) double {mustBePositive} = 48000;

    options.options_mu_gp  = gp_mu_opts;
    options.options_cov_gp = gp_cov_opts;

    options.options_disp_gp = gp_disp_opts('disp_legend_samples', false, 'disp_legend_loc', 'southoutside', 'disp_position', [100, 100, 800, 800], 'disp_legend_num_cols', 2);
    options.max_odr_gp (1,1) double {mustBeNonnegative} = 3;
    options.sample_method_gp (1,:) char {mustBeMember(options.sample_method_gp, {'mvnrnd', 'mean'})} = 'mean';
    options.N_freq_val_gp (1,1) double {mustBePositive, mustBeInteger} = 32;

    options.N_freq_uni_fit  (1,1) double {mustBePositive, mustBeInteger} = 16;
    options.max_taps_g (1,1) double {mustBePositive, mustBeInteger}  = 16;    

    options.SH_fit_mode (1,:) char {mustBeMember(options.SH_fit_mode, {'svd_inv', 'svd_ls'} )} = 'svd_inv';
    options.SH_fit_svd_trunc_frac (1,1) double {mustBeNonnegative} = 0;

    options.enable_disp (1,1) logical = false;
end

[N_C, M] = size(C);
P_C = sqrt(N_C) - 1;
assert(P_C == floor(P_C), 'Invalid size C');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create evaluation grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_D = options.max_odr_gp;
P_E = P_C + P_D;    %Combined max-order
N_E = (P_E + 1)^2;  %Total number of bases and evaluation points

[theta, phi] = sh_fib(N_E);
H = sh_dec(C, theta, phi, is_real); %[N_E x M]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample T60(omega, theta, phi) from GP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq = logspace(log10(20), log10(options.Fs/2), options.N_freq_val_gp)';
omega = freq * 2 * pi;
num_evals = 1;

[log_T60, h_fig_gp] = gp_t60_sample(omega, theta, phi, num_evals, obs, ...
    'options_mu', options.options_mu_gp, 'options_cov', options.options_cov_gp, ...
    'sample_method', options.sample_method_gp, ...
    'enable_disp', options.enable_disp, 'options_disp', options.options_disp_gp);
%log_T60 [options.N_freq_val x N_E x num_evals]

% [log_T60, h_fig_gp] = gp_t60_sample(omega, pi/2 * ones(7, 1), deg2rad(linspace(0, 180, 7)'), num_evals, obs, ...
%     'options_mu', options.options_mu_gp, 'options_cov', options.options_cov_gp, ...
%     'sample_method', options.sample_method_gp, ...
%     'enable_disp', options.enable_disp, 'options_disp', options.options_disp_gp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Interpolate log_T60 over uniform frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq_uniform = linspace(0, options.Fs/2, options.N_freq_uni_fit)'; %DC to Nyquist

RT60_sec = exp(interp1(freq, log_T60, freq_uniform));       %[options.N_uniform_freq x N_E]
%Extrapolate to nearest neighbor
mask_min_outside = freq_uniform < min(freq);
mask_max_outside = freq_uniform > max(freq);
if any(mask_min_outside)
    RT60_sec(mask_min_outside, :) = exp(log_T60(1, :));
end
if any(mask_max_outside)
    RT60_sec(mask_max_outside, :) = exp(log_T60(end, :));
end
RT60_sec = real(RT60_sec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fit filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set targets
H_tgt_dB_DC_NQ  = -60 ./ (options.Fs * RT60_sec);               %Target magnitude dB
H_tgt_abs_DC_NQ = db2mag(H_tgt_dB_DC_NQ);                       %Target magnitude modulus  
X_mag = [H_tgt_abs_DC_NQ; flipud(H_tgt_abs_DC_NQ(2:end-1, :))];    

N = min(options.max_taps_g, 2 * (options.N_freq_uni_fit-1) ); %Num filter taps of g
g_mat = zeros([N, N_E]);
options_coneprog = optimoptions("coneprog", MaxIterations=300, ConstraintTolerance=1e-8, OptimalityTolerance=1e-6);
for n = 1:N_E
    g_mat(:, n) = ft_bnd_minphase(X_mag(:, n), N, 'options_coneprog', options_coneprog, ...
        'mode', 'least_squares', 'enable_disp', false);
    ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exponentiate filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = zeros(N_E, N * M);
for n = 1:N_E
    F(n, :) = ft_exp_conv_opt(H(n, :), g_mat(:, n).');
    n

    %plot_RIR(F(132, 1:8000)', plot_RIR_opts('Fs', options.Fs, 'clim', [-180, -80], 'win_size', 512, 'N_FFT', 512, 'colormap', 'hot'));
    %plot_RIR(F(32, 1:8000)', plot_RIR_opts('Fs', options.Fs, 'clim', [-180, -80], 'win_size', 512, 'N_FFT', 512, 'colormap', 'hot'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Transform back into SH domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y = sh_val(P_E, theta, phi, is_real);
% E = Y \ F;
if strcmp(options.SH_fit_mode, 'svd_inv')

    E = sh_fit_svd(F, theta, phi, P_E, is_real, options.SH_fit_svd_trunc_frac);

elseif strcmp(options.SH_fit_mode, 'svd_ls')

    E = sh_fit_svd(F, theta, phi, P_C, is_real, options.SH_fit_svd_trunc_frac);

else
    error('Unknown options.SH_fit_mode');
end

% dB_lim =  [-160, -80];
% win_size = 64;
% N_FFT = 512;

% for phi_disp = 0:60:180
%     f_phi = sh_dec(E, pi/2, deg2rad(phi_disp), is_real);
%     plot_RIR(f_phi(1:8000), plot_RIR_opts('Fs', options.Fs, 'clim', dB_lim, 'win_size', win_size, 'N_FFT', N_FFT, 'colormap', 'hot', 'disp_RIR', false));
% end

% for phi_disp = 0:60:180
%     f_phi = sh_dec(C, pi/2, deg2rad(phi_disp), is_real);
%     plot_RIR(f_phi(1:8000), plot_RIR_opts('Fs', options.Fs, 'clim', dB_lim + 20, 'win_size', win_size, 'N_FFT', N_FFT, 'colormap', 'hot', 'disp_RIR', false));
% end

% Find closest angle to (pi/2, pi)
%V = zeros(N_E, 3);
%[V(:, 1), V(:, 2), V(:,3)] = sph2cart(phi, pi/2 - theta, ones(size(theta)));
%[u(1), u(2), u(3)] = sph2cart(deg2rad(180), 0, 1);
%[~, idx] = max(u*V')

%f_u = sh_dec(E, theta(idx), phi(idx), is_real);
%plot_RIR(f_u(1:8000), plot_RIR_opts('Fs', options.Fs, 'clim', dB_lim, 'win_size', win_size, 'N_FFT', N_FFT, 'colormap', 'hot', 'disp_RIR', false));

%f_u_perturb = sh_dec(E, theta(idx), phi(idx) + deg2rad(1), is_real);
%plot_RIR(f_u_perturb(1:8000), plot_RIR_opts('Fs', options.Fs, 'clim', dB_lim, 'win_size', win_size, 'N_FFT', N_FFT, 'colormap', 'hot', 'disp_RIR', false));

;