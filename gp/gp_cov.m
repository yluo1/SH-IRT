function [K, dK_mat_list, dK_name_list] = gp_cov(X, Y, options)
%Compute GP prior covariance

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X:        [NX x 3]     Row-matrix [omega, theta, phi]
%Y:        [NY x 3]     Row-matrix [omega, theta, phi]
%options:  Struct, see gp_cov_opts.m 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%K:        [NX x NY] Covariance matrix

%dK_mat_list:  [1 x P] cell matrix where dK_mat_list{1, n} is partial derivative matrix
%dK_name_list: [1 x P] cell matrix where dK_name_list{1, n} is name of variable in gp_cov_opts.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate sample covariance matrix and partial derivatives

% rng(234);
% N_B = 16;
% N_E = 10;
% Fs = 48000;
% omega = logspace(log10(20), log10(Fs/2), N_B)' * 2 *pi;
% [theta, phi] = sh_fib(N_E);
% freq = max(1, omega / (2 * pi)); %[N_B x 1]
% theta = rand(N_E, 1) * pi;
% phi = rand(N_E, 1) * 2 * pi;
% freq_grid = repmat(freq, [1, N_E]);
% theta_grid  = repmat(theta', [N_B, 1]);
% phi_grid    = repmat(phi', [N_B, 1]);
% X = [2 * pi * freq_grid(:), theta_grid(:), phi_grid(:)]; %[N_B * N_E x 1]
% [K, dK_mat_list, dK_name_list] = gp_cov(X, X);

arguments
    X (:,3) double = [1 0 0];
    Y (:,3) double = [1 0 0];
    options = gp_cov_opts;
end

if strcmp(options.cov_func, 'cov_sqx_chw_ns')
       
    if nargout > 1
        dK_mat_list = cell(1, 3);
        [K, dK_mat_list{:}] = cov_sqx_chw_ns(X, Y, options.cov_sigma, options.cov_ell, options.cov_gamma); %[NX x NY] K

        dK_name_list = {'cov_sigma', 'cov_ell', 'cov_gamma'};

    else
        K = cov_sqx_chw_ns(X, Y, options.cov_sigma, options.cov_ell, options.cov_gamma); %[NX x NY] K
    end

elseif strcmp(options.cov_func, 'cov_sqx_chw_sqx')

    if nargout > 1
        dK_mat_list = cell(1, 3);
        [K, dK_mat_list{:}] = cov_sqx_chw_sqx(X, Y, options.cov_sigma, options.cov_ell_c, options.cov_ell_f); %[NX x NY] K

        dK_name_list = {'cov_sigma', 'cov_ell_c', 'cov_ell_f'};

    else
        K = cov_sqx_chw_sqx(X, Y, options.cov_sigma, options.options.cov_ell_c, options); %[NX x NY] K
    end


else
    error('Unsupported options.cov_func');
end