function [err, h_ref, h_fit] = plot_sh_fit_example(mode, err_name, options)
%Plot sh_fit_* modes

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%mode:          String, fitting method {'TSVD', 'TR', 'SqExp', 'Mat52', 'Mat32', 'Exp'}
%err_name:      String, error function {'norm', 'SHMSQ', 'NSHMSQ'}

%options:       struct

%options.max_odr:           Maximum SH order for sample data generation   
%options.max_odr_fit:       Maximum SH order for fitting

%options.N_pts:             Number of sample spherical coordinates
%options.sample_mode:       String, sampling method {'fib', 'rand'}
%                               'fib':  Fibonnaci spiral
%                               'rand': Uniform random on sphere
%options.noise_std:         Noise standard deviation added to samples
%options.rseed:             Random seed

%options.svd_mode:          String, SVD regularization mode {'max', 'totalvar', 'bottom', 'picard'} (see sh_fit_svd.m)
%options.svd_trunc_frac:    Fraction of largest singular values to truncate in 'TSVD' mode

%options.tr_mode:           String, Tikhonov regularization mode {'identity', 'quad', 'quadlin', 'picard'} (see sh_fit_tr.m)
%options.tr_trunc_frac:     Regularization parameter in 'tr' mode

%options.rbf_max_iter:      Maximum number of hyper parameter optimization iterations for RBF modes
%options.rbf_objective:     String, hyper parameter optimization objective  {'NLMH', 'MSE'} (See sh_fit_rbf.m)
%options.rbf_ell:           Scalar, spatial bandwidth hyper-parameter       
%options.rbf_lambda:        Scalar, noise variance

%options.enable_disp:       Logical, if true, plot reference and target fields
%options.disp_dB_lim:       [1 x 2] Display dB range [min, max]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%err:       Error
%h_ref:     Figure handle to reference
%h_fit:     Figure handle to fit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Evaluate and plot sample fits

%[err_svd]      = plot_sh_fit_example('TSVD', 'NSHMSQ', 'svd_trunc_frac', 0.1)
%[err_SqExp0]   = plot_sh_fit_example('SqExp', 'NSHMSQ', 'rbf_max_iter', 0)

%[err_SqExp]    = plot_sh_fit_example('SqExp', 'NSHMSQ')
%[err_Exp]      = plot_sh_fit_example('Exp', 'NSHMSQ')

%varargin       = {'max_odr', 5, 'max_odr_fit', 10, 'N_pts', 100};
%[err_svd]      = plot_sh_fit_example('TSVD', 'NSHMSQ', varargin{:}, 'svd_mode', 'picard', 'svd_trunc_frac', 0.5)

%[err_tr]       = plot_sh_fit_example('TRegu', 'NSHMSQ', varargin{:}, 'tr_mode', 'picard', 'tr_lambda', 0.5)

%[err_SqExp]    = plot_sh_fit_example('SqExp', 'NSHMSQ', varargin{:})
%[err_Mat52]    = plot_sh_fit_example('Mat52', 'NSHMSQ', varargin{:})
%[err_Mat32]    = plot_sh_fit_example('Mat32', 'NSHMSQ', varargin{:})
%[err_Exp]      = plot_sh_fit_example('Exp', 'NSHMSQ', varargin{:})

arguments
    mode (1,:) char {mustBeMember(mode, {'TSVD', 'TRegu', 'SqExp', 'Mat52', 'Mat32', 'Exp'} )} = 'TSVD';
    err_name (1,:) char {mustBeMember(err_name, {'norm', 'SHMSQ', 'NSHMSQ' })} = 'NSHMSQ';

    options.max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 3;
    options.max_odr_fit (1,1) double {mustBeNonnegative, mustBeInteger} = 6;

    options.N_pts (1,1) double {mustBePositive, mustBeInteger} = 20;
    options.sample_mode (1,:) char {mustBeMember(options.sample_mode, {'fib', 'rand'} )} = 'rand';
    options.noise_std (1,1) double {mustBeNonnegative} = 1e-1;
    options.rseed (1,1) double {mustBeNonnegative, mustBeInteger} = 4141;

    options.svd_mode (1,:) char {mustBeMember(options.svd_mode, {'max', 'totalvar', 'bottom', 'picard'} )} = 'max';
    options.svd_trunc_frac (1,1) double {mustBeNonnegative} = 0;

    options.tr_mode (1,:) char {mustBeMember(options.tr_mode, {'identity', 'quad', 'quadlin', 'picard'} )} = 'identity';
    options.tr_lambda (1,1) double {mustBeNonnegative} = 0;

    options.rbf_max_iter (1,1) double {mustBeNonnegative, mustBeInteger} = 100;
    options.rbf_objective (1,:) char {mustBeMember(options.rbf_objective, {'NLMH', 'MSE'})} = 'NLMH';
    options.rbf_ell (1,1) double {mustBeNonnegative} = 0.5;
    options.rbf_lambda (1,1) double {mustBeNonnegative} = 0;

    options.enable_disp (1,1) logical = true;
    options.disp_dB_lim (1,2) double =  [-40, 20];
end

rng(options.rseed);
max_odr = options.max_odr;
is_real = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate random field and sample observations at spherical coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = sh_rand(max_odr, 1, is_real);
C = sh_resize(C, options.max_odr_fit);

% Sample spherical coordinates
if strcmp(options.sample_mode, 'fib')
    [theta, phi] = sh_grd_fib(options.N_pts);
elseif strcmp(options.sample_mode, 'rand')
    [theta, phi] = sh_grd_rand(options.N_pts);
else
    error('Unsupported options.sample_mode');
end

X = sh_dec(C, theta, phi, is_real);
X = X + (randn(size(X)) + randn(size(X)) * 1i) * options.noise_std; % Add noise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode, 'TSVD')

    C_fit = sh_fit_svd(X, theta, phi, options.max_odr_fit, is_real, options.svd_trunc_frac, options.svd_mode); 

elseif strcmp(mode, 'TRegu')

    C_fit = sh_fit_tr(X, theta, phi, options.max_odr_fit, is_real, options.tr_lambda, options.tr_mode);     

elseif any(strcmp(mode, {'SqExp', 'Mat52', 'Mat32', 'Exp'}))

    varargin = {'max_iter', options.rbf_max_iter, 'ub', [inf, 0], 'objective', options.rbf_objective};
    C_fit = sh_fit_rbf(X, theta, phi, options.max_odr_fit, is_real, mode, options.rbf_ell, options.rbf_lambda, varargin{:}); 

else
    error('Unsupported mode');
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(err_name, 'norm')

    err = norm(C - C_fit);

elseif strcmp(err_name, 'SHMSQ')

    err = err_SHMSQ(C, C_fit);

elseif strcmp(err_name, 'NSHMSQ')
    
    err = err_NSHMSQ(C, C_fit);

else
    error('Unsupported err_name');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_ref = [];
h_fit = [];
if options.enable_disp
    h_ref = sh_plt(C, 'mercator', is_real, 'dB_lim', options.disp_dB_lim, 'title_name', 'Reference', 'disp_theta_phi', [theta, phi]);
    h_fit = sh_plt(C_fit, 'mercator', is_real, 'dB_lim', options.disp_dB_lim, 'title_name', mode, 'disp_theta_phi', [theta, phi]); 

    h_ref = h_ref{1};
    h_fit = h_fit{1};
end
