function [C, err] = sh_fit_rbf(X, theta, phi, max_odr, is_real, mode, ell, lambda, options)
%Radial basis function kernel regression
%C = (K + lambda * I)^(-1) * (X - mean(X, 1));

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X:             [N x M] N measurements of M functions

%theta:         [N x 1] Co-latitude [0, pi]
%phi:           [N x 1] Azimuth [0, 2 * pi)

%max_odr:       Max SH order 

%is_real:       Logical, if true, evaluate real SH

%mode:          String, kernel function {'SqExp', 'Mat52', 'Mat32', 'Exp'}
%ell:           Spatial bandwidth
%lambda:        Noise variance

%options:               struct
%options.max_iter:      Maximum iterations for fitting parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:                 [(max_odr + 1)^2 x M] SH coefficients
%err:               Scalar, norm(Y*C - X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Fit to random field over random sampled points over sphere

% rng(4141);
% max_odr = 3;
% is_real = false;
% N_pts = (max_odr + 1)^2;
% C = sh_rand(max_odr, 1, is_real);
% %[theta, phi] = sh_fib(N_pts);
% [theta, phi] = sh_rand_unis(N_pts);
% X = sh_dec(C, theta, phi, is_real);
% X = X + (randn(size(X)) + randn(size(X)) * 1i) * 1e-1; % Add noise


% max_odr_fit = 6;
% ell = 0.5;
% lambda = 0;
% ub = [inf, 0.1];
% max_iter = 0;
% dB_lim = [-40, 20];
% C = sh_resize(C, max_odr_fit);

% sh_plt(C, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Reference', 'disp_theta_phi', [theta, phi]);
%
% [C_SqExp] = sh_fit_rbf(X, theta, phi, max_odr_fit, is_real, 'SqExp', ell, lambda, 'max_iter', max_iter, 'ub', ub); 
% err_SqExp = norm(C - C_SqExp)
% sh_plt(C_SqExp, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'SqExp', 'disp_theta_phi', [theta, phi]); 
% 
% [C_Mat52] = sh_fit_rbf(X, theta, phi, max_odr_fit, is_real, 'Mat52', ell, lambda, 'max_iter', max_iter, 'ub', ub); 
% err_Mat52 = norm(C - C_Mat52)
% sh_plt(C_Mat52, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Mat32', 'disp_theta_phi', [theta, phi]);
% 
% [C_Exp] = sh_fit_rbf(X, theta, phi, max_odr_fit, is_real, 'Exp', ell, lambda, 'max_iter', max_iter, 'ub', ub); 
% err_Exp = norm(C - C_Exp)
% sh_plt(C_Exp, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Exp', 'disp_theta_phi', [theta, phi]);

% [C_svd] = sh_fit_svd(X, theta, phi, max_odr_fit, is_real, 0); 
% err_svd = norm(C - C_svd)
% sh_plt(C_svd, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'SVD', 'disp_theta_phi', [theta, phi]);


arguments
    X (:,:) double {coder.mustBeComplex} = complex(0);

    theta (:,1) double = [0];
    phi   (:,1) double = [0];

    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;
    is_real (1,1) logical = false;

    mode (1,:) char {mustBeMember(mode, {'SqExp', 'Mat52', 'Mat32', 'Exp'})}  = 'SqExp';
    ell (1,1) double {mustBePositive} = 1;
    lambda (1,1) double {mustBeNonnegative} = 0;

    options.max_iter (1,1) double {mustBeNonnegative, mustBeInteger} = 100;
    options.lb (1,2) double {mustBeNonnegative} = [0, 0];
    options.ub (1,2) double {mustBeNonnegative} = [inf, inf];
    
end

[N, M] = size(X);

% Compute Chordal distance
VX = zeros([N, 3]);
[VX(:,1), VX(:,2), VX(:,3)] = sph2cart( phi, pi/2 - theta, ones(N, 1) );
D = pdist2(VX, VX); % Euclidean

X_mean = mean(X, 1);
X_centered = bsxfun(@minus, X, X_mean); % Center observations about mean

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximize log marginal likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.max_iter > 0
    options_fmincon = optimoptions("fmincon", SpecifyObjectiveGradient=false, Display="iter", checkGradients=false, ...
        ScaleProblem=true, ...
        FunctionTolerance=1e-8, ConstraintTolerance=1e-8, OptimalityTolerance=1e-10, StepTolerance=1e-8, ...
        MaxIterations=options.max_iter, MaxFunctionEvaluations=10000);
    
    param_list_0 = [ell, lambda];
    [param_list, fval, exitflag] = fmincon(@(x) neg_log_marginal_likelihood(x, mode, D, X_centered), ... 
        param_list_0, [], [], [], [], options.lb, options.ub, [], options_fmincon);
    
    ell = param_list(1);
    lambda = param_list(2);

    ell
    lambda
end

K = rbf_val(mode, D, ell);
K_inv_X = (K + lambda * eye(N)) \ X_centered;

C = sh_enc_rbf(mode, max_odr, theta, phi, ell, is_real) * K_inv_X;
C(1, :) = C(1, :) + X_mean * 2 * sqrt(pi); % Add mean

% Compute error
if nargout > 1
    Y   = sh_val(max_odr, theta, phi, is_real);
    err = norm(Y * C - X);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute log-marginal likelihood

%Input
%param_list:        [1 x 2] [ell, lambda]
%mode:              String, kernel function {'SqExp', 'Mat52', 'Mat32', 'Exp'}

%D:                 [N x N] Distance matrix
%X_centered:        [N x 1] Centered measurements

%Output
%fval:              negative log marginal likelihood

function fval = neg_log_marginal_likelihood(param_list, mode, D, X_centered)

N = numel(X_centered);

ell = param_list(1);
lambda = param_list(2);

K = rbf_val(mode, D, ell) + lambda * eye(N);

K_inv_X = K \ X_centered;

fval = 1/2 * ( trace(real( X_centered' * K_inv_X )) + sum(log( eig(K) )) + N * log(2*pi) );
