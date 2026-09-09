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

arguments
    omega (:,:) double {mustBeNonnegative} = pi/2;
    options = gp_mu_opts;
end

%Compute prior mean
if strcmp(options.mu_func, 'power')

    mu = mu_pow(omega, options.mu_alpha, options.mu_beta); %[N x M]

elseif strcmp(options.mu_func, 'LPF')

    mu = mu_lpf(omega, options.mu_alpha, options.mu_beta, options.mu_fc * 2 * pi); %[N x M]

else
    error('Unsupported options.mu_func');
end

