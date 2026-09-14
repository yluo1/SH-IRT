function mu = gp_mu(omega, options)
%Compute GP prior mean

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%omega:     [N x M] Angular frequency
%options:   Struct, see gp_mu_opts.m 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%mu:        [N x M]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Evaluate and plot prior mean function 

% omega = 2 * pi * logspace(log10(20), log10(24000), 256)'; % Angular frequencies

% mu_power = gp_mu(omega, gp_mu_opts('mu_func', 'power', 'mu_alpha', 1, 'mu_beta', 0.5, 'enable_disp', true));
% mu_lpf   = gp_mu(omega, gp_mu_opts('mu_func', 'LPF',   'mu_alpha', 1, 'mu_beta', 0.5, 'mu_fc', 500, 'enable_disp', true));

arguments
    omega (:,:) double {mustBeNonnegative} = pi/2;
    options = gp_mu_opts;
end

%Compute prior mean
if strcmp(options.mu_func, 'power')

    mu = mu_pow(omega, options.mu_alpha, options.mu_beta); %[N x M]

    if options.enable_disp  && coder.target('MATLAB')
        gp_plt('mu_pow', 'mu_pow_alpha_list', options.mu_alpha, 'mu_pow_beta_list', options.mu_beta, 'mu_pow_fig_position', [100, 100, 900, 450]);
    end
     
elseif strcmp(options.mu_func, 'LPF')

    mu = mu_lpf(omega, options.mu_alpha, options.mu_beta, options.mu_fc * 2 * pi); %[N x M]

    if options.enable_disp  && coder.target('MATLAB')
        gp_plt('mu_lpf', 'mu_lpf_alpha_list', options.mu_alpha, 'mu_lpf_beta_list', options.mu_beta, 'mu_lpf_omega_fc_list', options.mu_fc * 2 * pi, 'mu_lpf_fig_position', [100, 100, 900, 450]);
    end

else
    error('Unsupported options.mu_func');
end
