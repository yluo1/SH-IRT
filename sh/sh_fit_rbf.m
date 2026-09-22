function [C, err] = sh_fit_rbf(X, theta, phi, max_odr, is_real, mode, ell, lambda, options)
%Radial basis function (RBF) kernel expansion:
%C = K_exp *  (K + lambda * I)^(-1) * (X - mean(X, 1));

%where K_exp is the expansion of RBF kernels along SH bases [(max_odr + 1)^2 x N],
%K is the covariance matrix of RBFs [N x N]
%lambda is the noise variance term for diagonal loading

%Paper reference: 
%Luo, Y., 2021. Spherical harmonic covariance and magnitude function encodings for beamformer design.
%EURASIP Journal on Audio, Speech, and Music Processing, 2021(1), p.41.

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X:             [N x M] N measurements of M functions

%theta:         [N x 1] Co-latitude [0, pi]
%phi:           [N x 1] Azimuth [0, 2 * pi)

%max_odr:       Max SH order 

%is_real:       Logical, if true, evaluate real SH

%mode:          String, kernel function {'SqExp', 'SqExpNorm', 'Mat52', 'Mat32', 'Exp', 'ExpNorm'}
%ell:           Spatial bandwidth hyperparameter of covariance function
%lambda:        Noise variance

%options:               struct
%options.max_iter:      Maximum iterations for fitting parameters
%options.lb:            [1 x 2]     Lower bound [ell_min, lambda_min]
%options.ub:            [1 x 2]     Upper bound [ell_max, lambda_max]
%options.objective:     String, objective function {'NLMH', 'MSE'}
%                           'NLMH':     Negative log marginal likelihood
%                           'MSE':      Mean squared error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:                 [(max_odr + 1)^2 x M] SH coefficients
%err:               Scalar, norm(Y*C - X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: See plot_sh_fit_example.m

arguments
    X (:,:) double {coder.mustBeComplex} = complex(0);

    theta (:,1) double = [0];
    phi   (:,1) double = [0];

    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;
    is_real (1,1) logical = false;

    mode (1,:) char {mustBeMember(mode, {'SqExp', 'SqExpNorm', 'Mat52', 'Mat32', 'Exp', 'ExpNorm'})}  = 'SqExp';
    ell (1,1) double {mustBePositive} = 1;
    lambda (1,1) double {mustBeNonnegative} = 0;

    options.max_iter (1,1) double {mustBeNonnegative, mustBeInteger} = 100;
    options.lb (1,2) double {mustBeNonnegative} = [0, 0];
    options.ub (1,2) double {mustBeNonnegative} = [inf, inf];
    options.objective (1,:) char {mustBeMember(options.objective, {'NLMH', 'MSE'})} = 'NLMH';
    
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

    param_list_0 = [ell, lambda]';
    if strcmp(options.objective, 'NLMH')

        options_fmincon = optimoptions("fmincon", SpecifyObjectiveGradient=true, Display="final", checkGradients=false, ...
            ScaleProblem=true, ...
            FunctionTolerance=1e-8, ConstraintTolerance=1e-8, OptimalityTolerance=1e-10, StepTolerance=1e-8, ...
            MaxIterations=options.max_iter, MaxFunctionEvaluations=10000);

        func_obj = @(x) neg_log_marginal_likelihood(x, mode, D, X_centered);

    elseif strcmp(options.objective, 'MSE')

        options_fmincon = optimoptions("fmincon", SpecifyObjectiveGradient=false, Display="iter", checkGradients=false, ...
            ScaleProblem=true, ...
            FunctionTolerance=1e-8, ConstraintTolerance=1e-8, OptimalityTolerance=1e-10, StepTolerance=1e-8, ...
            MaxIterations=options.max_iter, MaxFunctionEvaluations=10000);

        func_obj = @(x) mean_squared_error(x, theta, phi, max_odr, is_real, mode, ...
            D, X_centered, sh_val(max_odr, theta, phi, is_real));

    else
        error('Unsupported options.objective')
    end
    [param_list, fval, exitflag] = fmincon(func_obj, ... 
        param_list_0, [], [], [], [], options.lb, options.ub, [], options_fmincon);
    
    ell     = param_list(1)
    lambda  = param_list(2)

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
%Compute negative log-marginal likelihood

%Input
%param_list:        [2 x 1] [ell, lambda]
%mode:              String, kernel function {'SqExp', 'Mat52', 'Mat32', 'Exp'}

%D:                 [N x N] Distance matrix
%X_centered:        [N x M] Centered measurements

%Output
%fval:              negative log marginal likelihood
%grad:              [2 x 1] df / dell

function [fval, grad] = neg_log_marginal_likelihood(param_list, mode, D, X_centered)

N = numel(X_centered);

ell = param_list(1);
lambda = param_list(2);

if nargout > 1
    [K, dK_dell] = rbf_val(mode, D, ell);
    dK_dlambda   = eye(N);
else
    K = rbf_val(mode, D, ell);
end

K = K + lambda * eye(N);

K_inv_X = K \ X_centered; %[N x M]

fval = 1/2 * ( trace(real( X_centered' * K_inv_X )) + sum(log( eig(K) )) + N * log(2*pi) );
fval = real(fval);

if nargout > 1
    grad    = zeros(2, 1);
    grad(1) = 1/2 * (trace( (K \ dK_dell) )    - K_inv_X' * dK_dell * K_inv_X );
    grad(2) = 1/2 * (trace( (K \ dK_dlambda) ) - K_inv_X' * dK_dlambda * K_inv_X );
    grad = real(grad);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute mean squared error

%Input
%param_list:        [1 x 2] [ell, lambda]
%mode:              String, kernel function {'SqExp', 'Mat52', 'Mat32', 'Exp'}

%D:                 [N x N] Distance matrix
%X_centered:        [N x M] Centered measurements

%X_mean:            [1 x M] Means per function
%Y:                 [N x (max_ord + 1)^2] SH bases evaluated over theta, phi

%Output
%fval:              negative log marginal likelihood

function fval = mean_squared_error(param_list, theta, phi, max_odr, is_real, mode, D, X_centered, Y)

N = numel(X_centered);

ell = param_list(1);
lambda = param_list(2);

K = rbf_val(mode, D, ell) + lambda * eye(N);

K_inv_X = K \ X_centered;

C = sh_enc_rbf(mode, max_odr, theta, phi, ell, is_real) * K_inv_X;

% Compute error
fval = mean(mean(abs(Y * C - X_centered ).^2));
