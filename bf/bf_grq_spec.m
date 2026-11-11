function [W] = bf_grq_spec(A, R, tau, mode, options)
% Generalized Rayleigh quotient beamformer designs for spherical harmonic expansion based 
% generalized directivity index GDI(A, R, W) = GRQ(A, R, W) =  w'*A*w / w'*R*w = tau,
% that minimize weighted spectral averages of the gengeralized eigenvalues

% tau truncates to the range [min(w'*A*w / w'*R*w), max(w'*A*w / w'*R*w)]

% mode minimizes different objective weights wsq_n = |w_n|^2 s.t.  subject to spectral GDI constraints:
% \sum_n wsq_n = 1,   wsq_n  >= 0,    Barycentric constraints
% \sum_n wsq_n * lambda_n = tau,      GDI constraints
% Rayleigh quotient solution: Weighted eigenvectors y = V * sqrt(wsq)

% mode = 'max':         Minimax squared barycenter weights (infinity norm) of GRQ generalized eigenvalues
% min_wsq max_n wsq_n,       

% mode = 'sq'           Minimize quintic barycenter weights (spectral norm) of GRQ generalized eigenvalues
% min_wsq \sum_n wsq_n^2,

% mode = 'wsq'          Minimize squared spectral weighted barycenter weights (spectral norm) of GRQ generalized eigenvalues
% min_wsq \sum_n wsq_n^2 * lambda_n^2, 

% Optional constraints:
% Null space:    U'*W_min_max = 0

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%A:             [N x N]     Accept covariance matrix
%R:             [N x N]     Reject positive definite covariance matrix
%tau:           Scalar, target GDI 
%mode:          String, averaging method {'minimax'}

%options:       struct

%options.U:     Empty or [N x N_U]  Null space vectors
%               Disable constraints if empty

%options.check_constraints: Logical, if true check for constraints satisfication of solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%W:     [N x 1] Beamformer weights (conjugated), scaled to have unit norm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compute generalized Rayleigh quotient for varied modes

% rng(5464);
% N = 4;
% A = randn(N) + randn(N) * 1i;     A = A * A';
% R = randn(N) + randn(N) * 1i;     R = R * R';
% tau = 10;

% [W_max] = bf_grq_spec(A, R, tau, 'max', 'check_constraints', true)
% [W_sq] = bf_grq_spec(A, R, tau, 'sq', 'check_constraints', true)
% [W_wsq] = bf_grq_spec(A, R, tau, 'wsq', 'check_constraints', true)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compute generalized Rayleigh quotient for varied modes, and homogeneous constraints

% rng(5464);
% N = 8;
% A = randn(N) + randn(N) * 1i;     A = A * A';
% R = randn(N) + randn(N) * 1i;     R = R * R';
% U = randn(N,1) + randn(N,1) * 1i;
% tau = 5;

% [W_max] = bf_grq_spec(A, R, tau, 'max', 'U', U, 'check_constraints', true)
% [W_sq] = bf_grq_spec(A, R, tau, 'sq', 'U', U, 'check_constraints', true)
% [W_wsq] = bf_grq_spec(A, R, tau, 'wsq', 'U', U, 'check_constraints', true)


arguments
    A (:,:) double {coder.mustBeComplex} = complex(1);
    R (:,:) double {coder.mustBeComplex} = complex(1);
    tau (1,1) double = inf;
    mode (1,:) char {mustBeMember(mode, {'max', 'sq', 'wsq'} )} = 'max';

    options.U (:,:) double {coder.mustBeComplex} = [];
    options.check_constraints (1,1) logical = false;
end

assert(isequal(size(A), size(R)), 'A, R size mismatch');
N = size(A, 1);
assert(size(A, 2) == N, 'A, R must be square matrices');
N_U = size(options.U, 2);

[L, flag] = chol(R, 'lower');
assert(flag == 0, 'R is not positive definite');

Q = (L \ A) / L';          % Transform numerator: inv(L) * A * inv(L')

if ~isempty(options.U) % Transform variables via projection: P = I - U * inv(U'*U) * U', where W = P*y, U'*W = U'*P*y = 0

    assert(size(options.U, 1) == N, 'U and A size mismatch');
    Y = L \ options.U; % [N x N_U] 
    P = eye(N) - Y * ((Y'*Y) \ (Y') );
    Q = P * Q * P'; 

end

[V, D] = eig(Q);
lambda = diag(D);
[~, idx] = sort(lambda, 'descend');  idx = idx(1:(N-N_U));  % Sort and truncate eigenvalues and eigenvectors    
lambda = real(lambda(idx));    
V = V(:, idx);
tau = max(min(tau, lambda(1)), lambda(end)); % Truncate tau to GRQ range

if strcmp(mode, 'max')

    w_sq = minimax_spec(lambda, tau);
    W = L' \ (V * sqrt(w_sq));

elseif strcmp(mode, 'sq')

    w_sq = min_spec_sq(lambda, tau);
    W = L' \ (V * sqrt(w_sq));

elseif strcmp(mode, 'wsq')

    w_sq = min_spec_wsq(lambda, tau);
    W = L' \ (V * sqrt(w_sq));

else
    error('Unknown mode');
end

% Normalize
W = bsxfun(@rdivide, W, sqrt(sum(W .* conj(W), 1)));

% Check constraints of solution
if options.check_constraints
    tol = 1e-6;
    pass_tau = abs(real( W' * A * W / ( W' * R * W ) ) - tau) <= tol

    if ~isempty(options.U)
        pass_U = abs(options.U' * W) <= tol
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w_sq = minimax_spec(lambda, tau, options)
% Minimax squared barycenter weights infinity norm) of GRQ generalized eigenvalues
% minimize_wsq max_wsq wsq,     ( wsq_n = |w_n|^2 ) s.t.   
% \sum_n wsq_n = 1,   wsq_n  >= 0,    Barycentric constraints
% \sum_n wsq_n * lambda_n = tau,      GDI constraints

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% lambda:   [N x 1]     Eigenvalues
% tau:      Scalar,     Target GDI in [min(lambda), max(lambda)] range
% options:  linear progrogram optimoptions 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% w_sq:     [N x 1]     Squared weights

arguments
    lambda (:,1) double = 1;
    tau (1,1) double = 1;
    
    options = optimoptions('linprog', 'Display', 'none', 'MaxIterations', 10000, ...
        'OptimalityTolerance', 1e-8, 'ConstraintTolerance', 1e-8);
end

assert(tau >= min(lambda) && tau <= max(lambda), 'tau must be in [min(lambda), max(lambda)]');

N = numel(lambda);

% Variable format [w_sq; epigraph t]
f = [zeros(N,1); 1];

A = [eye(N), -ones(N, 1)];
b = zeros(N, 1);

Aeq_bary = [ones(1,N), 0];
beq_bary = 1;

Aeq_grq = [abs(real(lambda(:)))', 0];
beq_grq = tau;

Aeq = [Aeq_bary; Aeq_grq];
beq = [beq_bary; beq_grq];

[w_aug, fval, exitflag] = linprog(f, A, b, Aeq, beq, zeros(N + 1, 1), ones(N + 1,1), options);
w_sq = w_aug(1:N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w_sq = min_spec_sq(lambda, tau, options)
% Minimize quintic barycenter weights (spectral norm) of GRQ generalized eigenvalues
% min_wsq \sum_n wsq_n^2,     (wsq_n = |w_n|^2) s.t.   
% \sum_n wsq_n = 1,   wsq_n  >= 0,    Barycentric constraints
% \sum_n wsq_n * lambda_n = tau,      GDI constraints
% Rayleigh quotient solution: Weighted eigenvectors y = V * sqrt(wsq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% lambda:   [N x 1]     Eigenvalues
% tau:      Scalar,     Target GDI in [min(lambda), max(lambda)] range
% options:  linear progrogram optimoptions 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% w_sq:     [N x 1]     Squared weights

arguments
    lambda (:,1) double = 1;
    tau (1,1) double = 1;
    
    options = optimoptions('quadprog', 'Display', 'none', 'MaxIterations', 10000, ...
        'OptimalityTolerance', 1e-8, 'FunctionTolerance', 1e-8, 'ConstraintTolerance', 1e-8);
end

assert(tau >= min(lambda) && tau <= max(lambda), 'tau must be in [min(lambda), max(lambda)]');

N = numel(lambda);

Aeq_bary = [ones(1,N)];
beq_bary = 1;

Aeq_grq = [abs(real(lambda(:)))'];
beq_grq = tau;

Aeq = [Aeq_bary; Aeq_grq];
beq = [beq_bary; beq_grq];

[w_sq, fval, exitflag] = quadprog(eye(N), zeros(N,1), [], [], Aeq, beq, zeros(N,1), ones(N,1), [], options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w_sq = min_spec_wsq(lambda, tau, options)
%  Minimize squared spectral weighted barycenter weights (spectral norm) of GRQ generalized eigenvalues
% min_wsq \sum_n wsq_n^2 * lambda_n^2,     (wsq_n = |w_n|^2) s.t.   
% \sum_n wsq_n = 1,   wsq_n  >= 0,    Barycentric constraints
% \sum_n wsq_n * lambda_n = tau,      GDI constraints
% Rayleigh quotient solution: Weighted eigenvectors y = V * sqrt(wsq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% lambda:   [N x 1]     Eigenvalues
% tau:      Scalar,     Target GDI in [min(lambda), max(lambda)] range
% options:  linear progrogram optimoptions 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% w_sq:     [N x 1]     Squared weights

arguments
    lambda (:,1) double = 1;
    tau (1,1) double = 1;
    
    options = optimoptions('quadprog', 'Display', 'none', 'MaxIterations', 10000, ...
        'OptimalityTolerance', 1e-8, 'FunctionTolerance', 1e-8, 'ConstraintTolerance', 1e-8);
end

assert(tau >= min(lambda) && tau <= max(lambda), 'tau must be in [min(lambda), max(lambda)]');

N = numel(lambda);

Aeq_bary = [ones(1,N)];
beq_bary = 1;

Aeq_grq = [abs(real(lambda(:)))'];
beq_grq = tau;

Aeq = [Aeq_bary; Aeq_grq];
beq = [beq_bary; beq_grq];

[w_sq, fval, exitflag] = quadprog(diag(lambda.^2), zeros(N,1), [], [], Aeq, beq, zeros(N,1), ones(N,1), [], options);


