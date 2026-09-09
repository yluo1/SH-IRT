function [E, h_fig_gp] = sh_exp_conv_gp(C, obs, options)
%Spherical harmonic exponentiating convolution of expansions C with
%direction-dependent g(theta, phi) sampled from Gaussian process T60(omega, theta, phi)

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:             [(P + 1)^2 x M]   SH coefficients (max order P of M number of functions)

%obs:                   Struct, observed T60(omega, theta, phi), see gp_obs_opts.m
%                       [] for none

%options.options_mu_gp:    Struct, prior mean options, see gp_mu_opts.m
%options.options_cov_gp:   Struct, prior covariance options, see gp_cov_opts.m

%options.max_odr_gp:       Max SH expansion order for functions drawn from GP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%E:  [(P + max_odr_gp + 1)^2 x M * max_taps_g]   SH coefficients
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
% E = sh_exp_conv_gp(C, obs, 'Fs', Fs, 'options_mu', options_mu, 'options_cov', options_cov, ...
%       'max_odr_gp', max_odr_gp, 'sample_method_gp', sample_method_gp, 'enable_disp', true);

% %Decode
% [theta, phi] = sh_plat('tetrahedron');
% f = real(sh_dec(E, theta, phi, is_real));
% plot_RIR(f(1, 1:Fs)');
% plot_RIR(f(2, 1:Fs)');

arguments
    C (:,:) double = 0;

    obs = [];

    options.Fs (1,1) double {mustBePositive} = 48000;
    options.is_real (1,1) logical = false;

    options.options_mu_gp  = gp_mu_opts;
    options.options_cov_gp = gp_cov_opts;

    options.options_disp_gp = gp_disp_opts('disp_legend_samples', false, 'disp_legend_loc', 'southoutside', 'disp_position', [100, 100, 800, 800], 'disp_legend_num_cols', 2);
    options.max_odr_gp (1,1) double {mustBeNonnegative} = 3;
    options.sample_method_gp (1,:) char {mustBeMember(options.sample_method_gp, {'mvnrnd', 'mean'})} = 'mvnrnd';
    options.N_freq_val_gp (1,1) double {mustBePositive, mustBeInteger} = 32;

    options.N_freq_uni_fit  (1,1) double {mustBePositive, mustBeInteger} = 16;
    options.max_taps_g (1,1) double {mustBePositive, mustBeInteger}  = 16;    
    
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
H = sh_dec(C, theta, phi, options.is_real); %[N_E x M]

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Transform back into SH domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = sh_val(P_E, theta, phi, options.is_real);
E = Y \ F;