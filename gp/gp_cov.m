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
%Sample usage: Generate and plot prior covariance matrix and partial derivatives

% N_B = 16;     % Number of log-uniform angular frequencies
% N_E = 20;     % Number of uniform spherical coordinates

% omega         = 2 * pi * logspace(log10(20), log10(24000), N_B)';
% [theta, phi]  = sh_fib(N_E);
% freq          = max(1, omega / (2 * pi)); %[N_B x 1]

% freq_grid     = repmat(freq, [1, N_E]);
% theta_grid    = repmat(theta', [N_B, 1]);
% phi_grid      = repmat(phi', [N_B, 1]);

% X = [2 * pi * freq_grid(:), theta_grid(:), phi_grid(:)]; %[N_B * N_E x 1] Input vector

% [K, dK_mat_list, dK_name_list] = gp_cov(X, X, gp_cov_opts('cov_func', 'cov_sqx_chw_ns', 'cov_sigma', 1, 'cov_ell', 343, 'cov_gamma', 1, 'enable_disp', true));
% [K, dK_mat_list, dK_name_list] = gp_cov(X, X, gp_cov_opts('cov_func', 'cov_sqx_chw_sqx', 'cov_sigma', 1, 'cov_ell_c', 1, 'cov_ell_f', 2, 'enable_disp', true));

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

    if options.enable_disp  && coder.target('MATLAB')
        gp_plt('cov_sqx_chw_ns', 'cov_sqx_chw_ns_sigma', options.cov_sigma, 'cov_sqx_chw_ns_ell_list', options.cov_ell, 'cov_sqx_chw_ns_gamma_list', options.cov_gamma, 'cov_sqx_chw_ns_fig_position', [100, 100, 1200, 900]);
    end

elseif strcmp(options.cov_func, 'cov_sqx_chw_sqx')

    if nargout > 1
        dK_mat_list = cell(1, 3);
        [K, dK_mat_list{:}] = cov_sqx_chw_sqx(X, Y, options.cov_sigma, options.cov_ell_c, options.cov_ell_f); %[NX x NY] K

        dK_name_list = {'cov_sigma', 'cov_ell_c', 'cov_ell_f'};

    else
        K = cov_sqx_chw_sqx(X, Y, options.cov_sigma, options.cov_ell_c, options.cov_ell_f); %[NX x NY] K
    end

    if options.enable_disp  && coder.target('MATLAB')
        gp_plt('cov_sqx_chw_sqx', 'cov_sqx_chw_sqx_sigma', options.cov_sigma, 'cov_sqx_chw_sqx_ell_c_list', options.cov_ell_c, 'cov_sqx_chw_sqx_ell_f_list', options.cov_ell_f, 'cov_sqx_chw_sqx_fig_position', [100, 100, 1200, 900]);
    end

else
    error('Unsupported options.cov_func');
end
