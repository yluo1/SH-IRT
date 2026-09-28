function [C, err, h_figs] = sh_fit_tr(X, theta, phi, max_odr, is_real, lambda, mode, options)
%Tikhonov regularization least squares: min_C ||Y(theta, phi) * C - X||^2 + lambda * C'*Q*C
%C = (Y'Y + lambda * Q)^(-1) * Y' * X

%mode = 'identity'  unity w.r.t. SH degree
%Q = I

%mode = 'quad'   diagonal regularization with quadratic weighting w.r.t. SH degree
%Q = diag(q),   q = [0, ..., l^2, ...] for SH degree l

%mode = 'quadlin'   diagonal regularization with quadratic + linear weighting w.r.t. SH degree
%Q = diag(q),   q = [0, ..., l*(l+1), ...] for SH degree l

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X:             [N x M] N measurements of M functions

%theta:         [N x 1] Co-latitude [0, pi]
%phi:           [N x 1] Azimuth [0, 2 * pi)

%max_odr:       Max SH order 
%is_real:       Logical, if true, evaluate real SH

%trunc_frac:    Fraction of largest singular values
%mode:          String, regularization method {'quad', 'quadlin'}

%options:               struct
%options.W:             [], [N x 1], or [N x N] weighting matrix for weighted least squares
%                       [] uniform weighting
%                       [N x 1] Positive weight vector
%                           min_C ||diag(W.^(1/2)) * (Y(theta, phi) * C - X)||^2 + lambda * C'*Q*C
%                       [N x N] Positive definite weight matrix
%                           min_C ||chol(W) * (Y(theta, phi) * C - X)||^2 + lambda * C'*Q*C

%options.enable_disp:   Logical, if true, plot fit and Picard plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:                 [(max_odr + 1)^2 x M] SH coefficients
%err:               Scalar, norm(Y*C - X);
%h_figs:            Handle to figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare Tikhonov regularization modes for random field + noise

% rng(441);
% max_odr = 6;
% is_real = false;
% N_C = (max_odr + 1)^2;
% N_pts = (max_odr + 1)^2;
% 
% %Reference field
% C_ref = sh_rand(max_odr, 1, is_real);
% 
% %Sample at spherical coordinates
% [theta, phi] = sh_grd_rand(N_pts);  %Random over sphere
% 
% X = sh_dec(C_ref, theta, phi, is_real);
% %Add noise
% X = X + (randn(size(X)) + randn(size(X)) * 1i) * 5e-2;

%lambda = 0.1;
%C_ls = sh_fit_tr(X, theta, phi, max_odr, is_real, 0);
%C_quad = sh_fit_tr(X, theta, phi, max_odr, is_real, lambda, 'quad');
%C_quadlin = sh_fit_tr(X, theta, phi, max_odr, is_real, lambda, 'quadlin');
%C_identity = sh_fit_tr(X, theta, phi, max_odr, is_real, lambda, 'identity');

%dB_lim = [-40, 20];
%sh_plt(C_ref, 'mercator', is_real, 'title_name', 'Reference', 'dB_lim', dB_lim);
%sh_plt(C_ls, 'mercator', is_real, 'title_name', 'Non-Regularized Least Squares', 'dB_lim', dB_lim);
%sh_plt(C_identity, 'mercator', is_real, 'title_name', 'TR Least Squares identity', 'dB_lim', dB_lim);
%sh_plt(C_quad, 'mercator', is_real, 'title_name', 'TR Least Squares quaddeg', 'dB_lim', dB_lim);
%sh_plt(C_quadlin, 'mercator', is_real, 'title_name', 'TR Least Squares quadlindeg', 'dB_lim', dB_lim);


arguments
    X (:,:) double {coder.mustBeComplex} = complex(0);

    theta (:,1) double = [0];
    phi   (:,1) double = [0];

    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;

    is_real (1,1) logical = false;
    
    lambda (1,1) double {mustBeNonnegative} = 0;

    mode (1,:) char {mustBeMember(mode, {'identity', 'quad', 'quadlin'})} = 'identity';

    options.W (:, :) double = [];
    options.enable_disp (1,1) logical = false;
end

% Compute SH bases
Y = sh_val(max_odr, theta, phi, is_real);  %[N x (max_odr + 1)^2]
[N, N_C] = size(Y);

% Apply optional weighting matrix
if ~isempty(options.W)
    if numel(options.W) == N

        assert(all(options.W(:) > 0), 'options.W must be positive');
        W_sqrt = sqrt(options.W(:));
        Y = bsxfun(@times, Y, W_sqrt);
        X = bsxfun(@times, X, W_sqrt);

    elseif isequal(size(options.W),  [N, N])        

        [W_sqrt, flag] = chol(options.W);
        assert(flag == 0, 'options.W is not positive definite');
        Y = W_sqrt * Y;
        X = W_sqrt * X;

    else
        error('options.W size unsupported');
    end
end

% Compute regularization matrix
if strcmp(mode, 'identity')

    Q = eye(N_C);
    
elseif strcmp(mode, 'quad')

    q = zeros(N_C, 1);
    for l = 0:max_odr
        idx = (l + 1)^2 - l;
        q(idx) = l^2;
    end
    Q = diag(q);
    
elseif strcmp(mode, 'quadlin')

    q = zeros(N_C, 1);
    for l = 0:max_odr
        idx = (l + 1)^2 - l;
        q(idx) = l * (l + 1);
    end
    Q = diag(q);

else
    error('Unknown mode');
end

% Compute least squares solution
C = (Y'*Y + lambda * Q) \ (Y' * X);

% Compute error
err = norm(Y * C - X);

h_figs = [];
if options.enable_disp && coder.target("MATLAB")
    % Plot fit
    h_figs = sh_plt(C, 'mercator', is_real, 'disp_theta_phi', [theta, phi], 'title_name', ['TR Fit ', mode, ' \lambda = ', num2str(lambda)]);

end