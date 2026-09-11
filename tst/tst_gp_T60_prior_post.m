function tst_gp_T60_prior_post(mode, options)
%Generate test examples for sampling T60 from GP prior, posterior evaluation grids

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%mode:     String, test mode {'prior'}

%options:           Struct
%options.Fs:        Sampling rate
%options.rseed:     Random seed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:
%tst_gp_T60_prior_post('prior');

arguments
    mode (1,:) char {mustBeMember(mode, {'prior'})} = 'prior';    
    
    options.Fs (1,1) double {mustBePositive} = 48000;
    options.rseed (1,1) double {mustBePositive, mustBeInteger} = 21136;
end

if strcmp(mode, 'prior')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Setup GP prior mean and covariance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    options_mu = gp_mu_opts('mu_func', 'power', 'mu_alpha', 1, 'mu_beta', 0.25, 'mu_fc', 8000);
    options_cov = gp_cov_opts('cov_sigma', sqrt(2)/2, 'cov_gamma', 2/3, 'cov_ell', 343 * 1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Single point spherical coordinate evaluation grid, sampled 4 times, GP prior
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Setup evaluation grid
    N_B = 128;
    N_E = 1;
    omega = logspace(log10(20), log10(options.Fs/2), N_B)' * 2 *pi;
    [theta, phi] = sh_fib(N_E);
    num_evals = 4; %Sample 4 functions
    
    %Sample function
    rng(options.rseed + 8);
    [log_T60, h_fig_prior] = gp_t60_sample(omega, theta, phi, num_evals, [], ...
            'options_mu', options_mu, 'options_cov', options_cov, ...
            'enable_disp', true, 'options_disp', gp_disp_opts('disp_ylim', [0, 1.75], 'disp_legend_num_cols', 1));

else

    error('Unknown mode');

end