function [W_min_max, tau_min_max] = bf_grq_ex(A, R, options)
% Generalized Rayleigh quotient extrema beamformer design:
% tau = w'*A*w / w'*R*w
% where R is positive definite [N x N] matrix

% Paper references:
% Gander, W., Golub, G. H., & Von Matt, U. (1989). A constrained eigenvalue problem. Linear Algebra and its applications, 114, 815-839.
% Golub, G. H. (1973). Some modified matrix eigenvalue problems. SIAM review, 15(2), 318-334.
% Yunshen Zhou, Zhaojun Bai, Ren-Cang Li. Linear Constrained Rayleigh Quotient Optimization: Theory and Algorithms. CSIAM Trans. Appl. Math., 2021, 2 (2) : 195-262 DOI:10.4208/csiam-am.2021.nla.01

% Optional constraints:
% Null space:    U'*W_min_max = 0

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%A:             [N x N]     Accept covariance matrix
%R:             [N x N]     Reject covariance matrix

%options:       struct
%options.mode:  String, compute method {'chol', 'gev'}
%                   'chol':     Cholesky decomposition of R
%                   'gev':      Generalized eigenvector A*V = R*V*D

%options.U:     Empty or [N x N_U]  Null space vectors
%               Disable constraints if empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%W_min_max:     [N x 2] Critical points [minimizer, maximizer] scaled to have unit norm
%tau_min_max:   [1 x 2] Extrema [min(w'*A*w / w'*R*w), max(w'*A*w / w'*R*w)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compute generalized Rayleigh quotient extrema

% rng(5464);
% N = 4;
% A = randn(N) + randn(N) * 1i;     A = A * A';
% R = randn(N) + randn(N) * 1i;     R = R * R';

% [W_min_max_chol, tau_min_max_chol] = bf_grq_ex(A, R, 'mode', 'chol')
% [W_min_max_gev, tau_min_max_gev] = bf_grq_ex(A, R, 'mode', 'gev')

% err = norm(bsxfun(@rdivide, W_min_max_chol, W_min_max_chol(2, :)) - bsxfun(@rdivide, W_min_max_gev, W_min_max_gev(2, :))) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample usage: GRQ extrema with null space constraints

% rng(5464);
% N = 4;
% A = randn(N) + randn(N) * 1i;     A = A * A';
% R = randn(N) + randn(N) * 1i;     R = R * R';
% U = randn(N,1) + randn(N,1) * 1i;

% [W_min_max_chol, tau_min_max_chol] = bf_grq_ex(A, R, 'mode', 'chol', 'U', U);
% U' * W_min_max_chol
% tau_min_max_chol 

% [W_min_max_gev, tau_min_max_gev] = bf_grq_ex(A, R, 'mode', 'gev', 'U', U);
% U' * W_min_max_gev
% tau_min_max_gev

% err = norm(bsxfun(@rdivide, W_min_max_chol, W_min_max_chol(2, :)) - bsxfun(@rdivide, W_min_max_gev, W_min_max_gev(2, :))) 

arguments
    A (:,:) double {coder.mustBeComplex} = complex(1);
    R (:,:) double {coder.mustBeComplex} = complex(1);

    options.mode (1,:) char {mustBeMember(options.mode, {'chol', 'gev'} )} = 'chol';
    options.U (:,:) double {coder.mustBeComplex} = [];
end

assert(isequal(size(A), size(R)), 'A, R size mismatch');
N = size(A, 1);
assert(size(A, 2) == N, 'A, R must be square matrices');

if strcmp(options.mode, 'chol')

    [L, flag] = chol(R, 'lower');
    assert(flag == 0, 'R is not positive definite');

    Q = (L \ A) / L';          % Transform numerator: inv(L) * A * inv(L')

    if ~isempty(options.U) % Transform variables via projection: P = I - U * inv(U'*U) * U', where W = P*y, U'*W = U'*P*y = 0

        assert(size(options.U, 1) == N, 'U and A size mismatch');
        N_U = size(options.U, 2);
        Y = L \ options.U; % [N x N_U] 
        P = eye(N) - Y * ((Y'*Y) \ (Y') );
        Q = P * Q * P'; 

        [y_min, tau_min] = eigs(Q, 1 + N_U, 'smallestreal');
        y_min = y_min(:, end);
        tau_min = tau_min(end);
    else
        [y_min, tau_min] = eigs(Q, 1, 'smallestreal');
    end

    [y_max, tau_max] = eigs(Q, 1, 'largestreal'); 

    tau_min_max = real([tau_min, tau_max]);
    W_min_max = L' \ [y_min, y_max];

elseif strcmp(options.mode, 'gev')
  
    if isempty(options.U)

        [V_min, d_min] = eigs(A, R, 1, 'smallestreal');
        [V_max, d_max] = eigs(A, R, 1, 'largestreal');
    
        assert( isfinite(d_min) && isfinite(d_max), 'R is not positive definite');
        tau_min_max = real([d_min, d_max]);
        W_min_max = [V_min, V_max];

    else % Transform variables via projection: P = I - U * inv(U'*U) * U', where W = P*y, U'*W = U'*P*y = 0

        assert(size(options.U, 1) == N, 'U and A size mismatch');

        U = options.U;
        P = eye(N) - U * ((U'*U) \ U');

        A_null = P' * A * P;
        R_null = P' * R * P;
                
        [V_min, d_min] = eigs(A_null, R_null, 1, 'smallestreal');
        [V_max, d_max] = eigs(A_null, R_null, 1, 'largestreal');
    
        assert( isfinite(d_min) && isfinite(d_max), 'R is not positive definite');
        tau_min_max = real([d_min, d_max]);
        W_min_max = P * [V_min, V_max];

    end

    
else
    error('Unknown options.mode');
end

% Normalize
W_min_max = bsxfun(@rdivide, W_min_max, sqrt(sum(W_min_max .* conj(W_min_max), 1)));

% Check
% real( W_min_max(:, 1)' * A * W_min_max(:, 1) / ( W_min_max(:, 1)' * R * W_min_max(:, 1) ) ) - tau_min_max(1)
% real( W_min_max(:, 2)' * A * W_min_max(:, 2) / ( W_min_max(:, 2)' * R * W_min_max(:, 2) ) ) - tau_min_max(2)
