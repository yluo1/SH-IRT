function [W_min_max, tau_min_max] = bf_grq_ex(A, R, options)
% Generalized Rayleigh quotient extrema beamformer design
% tau = w'*A*w / w'*R*w
% where R is positive definite [N x N] matrix

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%A:             [N x N]     Accept covariance matrix
%R:             [N x N]     Reject covariance matrix

%options:       struct
%options.mode:  String, compute method {'chol', 'gev'}
%                   'chol':     Cholesky decomposition of R
%                   'gev':      Generalized eigenvector A*V = R*V*D

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

% err = norm(bsxfun(@rdivide, W_min_max_chol, W_min_max_chol(1, :)) - bsxfun(@rdivide, W_min_max_gev, W_min_max_gev(1, :))) 

arguments
    A (:,:) double {coder.mustBeComplex} = complex(1);
    R (:,:) double {coder.mustBeComplex} = complex(1);
    options.mode (1,:) char {mustBeMember(options.mode, {'chol', 'gev'} )} = 'chol';
end

assert(isequal(size(A), size(R)), 'A, R size mismatch');

if strcmp(options.mode, 'chol')

    [L, flag] = chol(R, 'lower');
    assert(flag == 0, 'R is not positive definite');

    Q = (L \ A) / L';          %inv(L) * A * inv(L')
    [y_min, tau_min] = eigs(Q, 1, 'smallestreal');
    [y_max, tau_max] = eigs(Q, 1, 'largestreal'); 

    tau_min_max = real([tau_min, tau_max]);
    W_min_max = L' \ [y_min, y_max];

elseif strcmp(options.mode, 'gev')

    [V, D] = eig(A, R);
    d = diag(D);
    [d, idx] = sort(d, 'ascend');
    V = V(:, idx);

    assert(all(isfinite(d)), 'R is not positive definite');
    tau_min_max = real([d(1), d(end)]);
    W_min_max = [V(:, 1), V(:, end)];

else
    error('Unknown options.mode');
end

% Normalize
W_min_max = bsxfun(@rdivide, W_min_max, sqrt(sum(W_min_max .* conj(W_min_max), 1)));

% Check
% real( W_min_max(:, 1)' * A * W_min_max(:, 1) / ( W_min_max(:, 1)' * R * W_min_max(:, 1) ) ) - tau_min_max(1)
% real( W_min_max(:, 2)' * A * W_min_max(:, 2) / ( W_min_max(:, 2)' * R * W_min_max(:, 2) ) ) - tau_min_max(2)
