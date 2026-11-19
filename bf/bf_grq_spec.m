function [W] = bf_grq_spec(A, R, tau, mode, options)
% Generalized Rayleigh quotient based spectral minimizations with constant directivity beamformer designs
% GRQ(A, R, W) =  w'*A*w / w'*R*w = tau,

% W = V * sqrt(w_sq)    for eigenvectors V, and squared weights w_sq,

% where w_sq is minimizer to the following modes that 
% minimizes different objective weights wsq_n = |w_n|^2 s.t.  subject to spectral GRQ constraints:
% \sum_n wsq_n = 1,   wsq_n  >= 0,    Barycentric constraints
% \sum_n wsq_n * lambda_n = tau,      GRQ constraints
% Rayleigh quotient solution: Weighted eigenvectors y = V * sqrt(wsq)

% For non-empty channel options.wt, bias the spectral weights by the squared norm 
% of the channel weighted un-transformed eigenvectors (L' \ V)
% wev(n) = 1,                                                      empty wt
% wev(n) = 1 ./ ( (L' \ V(:,n))' * diag(wt.^2) * (L' \ V(:,n)) ),  non-empty wt

% mode = 'max':         Minimize maximum barycenter squared weights (infinity norm) of GRQ generalized eigenvalues
% min_wsq max_n wev_n * wsq_n, 

% mode = 'wabs'         Minimize spectral weighted barycenter weights of GRQ generalized eigenvalues
% min_wsq \sum_n wev_n * wsq_n * lambda_n,

% mode = 'iwabs'         Minimize inverse spectral weighted barycenter weights of GRQ generalized eigenvalues
% min_wsq \sum_n wev_n * wsq_n * lambda_n^(-1),

% mode = 'sq'           Minimize barycenter squared weights (spectral norm) of GRQ generalized eigenvalues
% min_wsq \sum_n wev_n * wsq_n^2,

% mode = 'wsq'          Minimize squared spectral weighted barycenter squared weights (spectral norm) of GRQ generalized eigenvalues
% min_wsq \sum_n wev_n * wsq_n^2 * lambda_n^2, 

% mode = 'iwsq'         Minimize inverse squared spectral weighted barycenter squared weights (spectral norm) of GRQ generalized eigenvalues
% min_wsq \sum_n wev_n * wsq_n^2 * lambda_n^(-2), 

% tau truncates to the range [min(w'*A*w / w'*R*w), max(w'*A*w / w'*R*w)]

% Optional constraints:
% Null space:    U'*W_min_max = 0

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%A:             [N x N]     Accept covariance matrix
%R:             [N x N]     Reject positive definite covariance matrix
%tau:           Scalar, target GRQ 
%mode:          String, averaging method {'max', 'wabs', 'iwabs', 'sq', 'wsq', 'iwsq'}

%options:       struct

%options.U:     Empty or [N x N_U]  Null space vectors
%               Disable constraints if empty

%options.check_constraints: Logical, if true check for constraints satisfication of solution

%options.wt:    Empty or [N x 1] channel weights

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%W:     [N x 1] Beamformer weights (conjugated), scaled to have unit norm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compute solutions for varied modes

% rng(5464);
% N = 4;
% A = randn(N) + randn(N) * 1i;     A = A * A';
% R = randn(N) + randn(N) * 1i;     R = R * R';
% tau = 10;

% [W_max] = bf_grq_spec(A, R, tau, 'max', 'check_constraints', true)
% [W_sq] = bf_grq_spec(A, R, tau, 'sq', 'check_constraints', true)
% [W_wsq] = bf_grq_spec(A, R, tau, 'wsq', 'check_constraints', true)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compute solutions for varied modes, and homogeneous constraints

% rng(5464);
% N = 8;
% A = randn(N) + randn(N) * 1i;     A = A * A';
% R = randn(N) + randn(N) * 1i;     R = R * R';
% U = randn(N,1) + randn(N,1) * 1i;
% tau = 5;

% [W_max] = bf_grq_spec(A, R, tau, 'max', 'U', U, 'check_constraints', true)
% [W_sq] = bf_grq_spec(A, R, tau, 'sq', 'U', U, 'check_constraints', true)
% [W_wsq] = bf_grq_spec(A, R, tau, 'wsq', 'U', U, 'check_constraints', true)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compute solutions for channel weightings

% rng(5464);
% N = 8;
% A = randn(N) + randn(N) * 1i;     A = A * A';
% R = randn(N) + randn(N) * 1i;     R = R * R';
% U = randn(N,1) + randn(N,1) * 1i;
% tau = 5;

% wt = rand(N,1);
% [W_max_wt]    = bf_grq_spec(A, R, tau, 'max', 'wt', wt, 'check_constraints', true)
% [W_max_U_wt]  = bf_grq_spec(A, R, tau, 'max', 'wt', wt, 'U', U,  'check_constraints', true)

arguments
    A (:,:) double {coder.mustBeComplex} = complex(1);
    R (:,:) double {coder.mustBeComplex} = complex(1);
    tau (1,1) double = inf;
    mode (1,:) char {mustBeMember(mode, {'max', 'wabs', 'iwabs', 'sq', 'wsq', 'iwsq'} )} = 'max';

    options.U (:,:) double {coder.mustBeComplex} = [];
    options.check_constraints (1,1) logical = false;

    options.wt (:,1) double  {mustBeNonnegative}  = []; 
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

if tau >= lambda(1)  % Solution is GRQ maximizer 

    W = L' \ V(:, 1);

elseif tau <= lambda(end) % Solution is GRQ minimizer

    W = L' \ V(:, end);

else

    % Compute channel-weighted eigenvector weights
    if isempty(options.wt)
        wev = ones(N-N_U,1);
    else
        T = L' \ V;
        wev = zeros(N-N_U,1);
        for n = 1:(N-N_U)
            wev(n) = real(T(:,n)' * diag(options.wt.^(2)) * T(:,n));
        end
        wev = 1 ./ wev;
        wev = wev / sum(wev);
    end

    if strcmp(mode, 'max')
    
        w_sq = minimax_spec(lambda, tau, wev);

    elseif strcmp(mode, 'wabs')
    
        w_sq = min_spec_wabs(lambda, tau, wev);

    elseif strcmp(mode, 'iwabs')
    
        w_sq = min_spec_iwabs(lambda, tau, wev);
    
    elseif strcmp(mode, 'sq')
    
        w_sq = min_spec_sq(lambda, tau, wev);
    
    elseif strcmp(mode, 'wsq')
    
        w_sq = min_spec_wsq(lambda, tau, wev);

    elseif strcmp(mode, 'iwsq')
    
        w_sq = min_spec_iwsq(lambda, tau, wev);

    else
        error('Unknown mode');
    end
    W = L' \ (V * sqrt(w_sq));

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
function w_sq = minimax_spec(lambda, tau, wev, options)
% Minimax squared barycenter weights infinity norm) of GRQ generalized eigenvalues
% minimize_wsq max_wsq wev .* wsq,     ( wsq_n = |w_n|^2 ) s.t.   
% \sum_n wsq_n = 1,   wsq_n  >= 0,    Barycentric constraints
% \sum_n wsq_n * lambda_n = tau,      GRQ constraints

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% lambda:   [N x 1]     Eigenvalues
% tau:      Scalar,     Target GRQ in [min(lambda), max(lambda)] range
% wev:      [N x 1]     Channel-eigenvalue weighting

% options:  linprog optimoptions 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% w_sq:     [N x 1]     Squared weights

arguments
    lambda (:,1) double = 1;
    tau (1,1) double = 1;
    wev (:,1) double = 1;
    
    options = optimoptions('linprog', 'Display', 'none', 'MaxIterations', 10000, ...
        'OptimalityTolerance', 1e-8, 'ConstraintTolerance', 1e-8);
end

assert(tau >= min(lambda) && tau <= max(lambda), 'tau must be in [min(lambda), max(lambda)]');
assert(numel(lambda) == numel(wev), 'wev size mismatch lambda');

N = numel(lambda);

% Variable format [w_sq; epigraph t]
f = [zeros(N,1); 1];

A = [diag(wev), -ones(N, 1)];
b = zeros(N, 1);

Aeq_bary = [ones(1,N), 0];
beq_bary = 1;

Aeq_grq = [abs(real(lambda(:)))', 0];
beq_grq = tau;

Aeq = [Aeq_bary; Aeq_grq];
beq = [beq_bary; beq_grq];

[w_aug, fval, exitflag] = linprog(f, A, b, Aeq, beq, zeros(N + 1, 1), [ones(N,1); inf], options);
w_sq = w_aug(1:N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w_sq = min_spec_wabs(lambda, tau, wev, options)
% Minimize quintic barycenter weights (spectral norm) of GRQ generalized eigenvalues
% min_wsq \sum_n wev_n * wsq_n * lambda_n,     (wsq_n = |w_n|^2) s.t.   
% \sum_n wsq_n = 1,   wsq_n  >= 0,    Barycentric constraints
% \sum_n wsq_n * lambda_n = tau,      GRQ constraints
% Rayleigh quotient solution: Weighted eigenvectors y = V * sqrt(wsq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% lambda:   [N x 1]     Eigenvalues
% tau:      Scalar,     Target GRQ in [min(lambda), max(lambda)] range
% wev:      [N x 1]     Channel-eigenvalue weighting

% options:  linprog optimoptions 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% w_sq:     [N x 1]     Squared weights

arguments
    lambda (:,1) double = 1;
    tau (1,1) double = 1;
    wev (:,1) double = 1;
    
    options = optimoptions('linprog', 'Display', 'none', 'MaxIterations', 10000, ...
        'OptimalityTolerance', 1e-8, 'ConstraintTolerance', 1e-8);
end

assert(tau >= min(lambda) && tau <= max(lambda), 'tau must be in [min(lambda), max(lambda)]');
assert(numel(lambda) == numel(wev), 'wev size mismatch lambda');

N = numel(lambda);

Aeq_bary = [ones(1,N)];
beq_bary = 1;

Aeq_grq = [abs(real(lambda(:)))'];
beq_grq = tau;

Aeq = [Aeq_bary; Aeq_grq];
beq = [beq_bary; beq_grq];

f = lambda;

[w_sq, fval, exitflag] = linprog(f, [], [], Aeq, beq, zeros(N,1), ones(N,1), options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w_sq = min_spec_iwabs(lambda, tau, wev, options)
% Minimize quintic barycenter weights (spectral norm) of GRQ generalized eigenvalues
% min_wsq \sum_n wev_n * wsq_n * lambda_n^(-1),     (wsq_n = |w_n|^2) s.t.   
% \sum_n wsq_n = 1,   wsq_n  >= 0,    Barycentric constraints
% \sum_n wsq_n * lambda_n = tau,      GRQ constraints
% Rayleigh quotient solution: Weighted eigenvectors y = V * sqrt(wsq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% lambda:   [N x 1]     Eigenvalues
% tau:      Scalar,     Target GRQ in [min(lambda), max(lambda)] range
% wev:      [N x 1]     Channel-eigenvalue weighting

% options:  linprog optimoptions 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% w_sq:     [N x 1]     Squared weights

arguments
    lambda (:,1) double = 1;
    tau (1,1) double = 1;
    wev (:,1) double = 1;
    
    options = optimoptions('linprog', 'Display', 'none', 'MaxIterations', 10000, ...
        'OptimalityTolerance', 1e-8, 'ConstraintTolerance', 1e-8);
end

assert(tau >= min(lambda) && tau <= max(lambda), 'tau must be in [min(lambda), max(lambda)]');
assert(numel(lambda) == numel(wev), 'wev size mismatch lambda');

N = numel(lambda);

Aeq_bary = [ones(1,N)];
beq_bary = 1;

Aeq_grq = [abs(real(lambda(:)))'];
beq_grq = tau;

Aeq = [Aeq_bary; Aeq_grq];
beq = [beq_bary; beq_grq];

f = 1 ./ lambda;

[w_sq, fval, exitflag] = linprog(f, [], [], Aeq, beq, zeros(N,1), ones(N,1), options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w_sq = min_spec_sq(lambda, tau, wev, options)
% Minimize quintic barycenter weights (spectral norm) of GRQ generalized eigenvalues
% min_wsq \sum_n wev_n * wsq_n^2,     (wsq_n = |w_n|^2) s.t.   
% \sum_n wsq_n = 1,   wsq_n  >= 0,    Barycentric constraints
% \sum_n wsq_n * lambda_n = tau,      GRQ constraints
% Rayleigh quotient solution: Weighted eigenvectors y = V * sqrt(wsq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% lambda:   [N x 1]     Eigenvalues
% tau:      Scalar,     Target GRQ in [min(lambda), max(lambda)] range
% wev:      [N x 1]     Channel-eigenvalue weighting

% options:  quadprog optimoptions 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% w_sq:     [N x 1]     Squared weights

arguments
    lambda (:,1) double = 1;
    tau (1,1) double = 1;
    wev (:,1) double = 1;
    
    options = optimoptions('quadprog', 'Display', 'none', 'MaxIterations', 10000, ...
        'OptimalityTolerance', 1e-8, 'FunctionTolerance', 1e-8, 'ConstraintTolerance', 1e-8);
end

assert(tau >= min(lambda) && tau <= max(lambda), 'tau must be in [min(lambda), max(lambda)]');
assert(numel(lambda) == numel(wev), 'wev size mismatch lambda');

N = numel(lambda);

Aeq_bary = [ones(1,N)];
beq_bary = 1;

Aeq_grq = [abs(real(lambda(:)))'];
beq_grq = tau;

Aeq = [Aeq_bary; Aeq_grq];
beq = [beq_bary; beq_grq];

[w_sq, fval, exitflag] = quadprog(diag(wev), zeros(N,1), [], [], Aeq, beq, zeros(N,1), ones(N,1), [], options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w_sq = min_spec_wsq(lambda, tau, wev, options)
%  Minimize squared spectral weighted barycenter weights (spectral norm) of GRQ generalized eigenvalues
% min_wsq \sum_n wev_n * wsq_n^2 * lambda_n^2,     (wsq_n = |w_n|^2) s.t.   
% \sum_n wsq_n = 1,   wsq_n  >= 0,    Barycentric constraints
% \sum_n wsq_n * lambda_n = tau,      GRQ constraints
% Rayleigh quotient solution: Weighted eigenvectors y = V * sqrt(wsq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% lambda:   [N x 1]     Eigenvalues
% tau:      Scalar,     Target GRQ in [min(lambda), max(lambda)] range
% wev:      [N x 1]     Channel-eigenvalue weighting

% options:  quadprog ptimoptions 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% w_sq:     [N x 1]     Squared weights

arguments
    lambda (:,1) double = 1;
    tau (1,1) double = 1;
    wev (:,1) double = 1;

    options = optimoptions('quadprog', 'Display', 'none', 'MaxIterations', 10000, ...
        'OptimalityTolerance', 1e-8, 'FunctionTolerance', 1e-8, 'ConstraintTolerance', 1e-8);
end

assert(tau >= min(lambda) && tau <= max(lambda), 'tau must be in [min(lambda), max(lambda)]');
assert(numel(lambda) == numel(wev), 'wev size mismatch lambda');

N = numel(lambda);

Aeq_bary = [ones(1,N)];
beq_bary = 1;

Aeq_grq = abs(real(lambda(:)))';
beq_grq = tau;

Aeq = [Aeq_bary; Aeq_grq];
beq = [beq_bary; beq_grq];

[w_sq, fval, exitflag] = quadprog(diag(wev .* (lambda.^2)), zeros(N,1), [], [], Aeq, beq, zeros(N,1), ones(N,1), [], options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w_sq = min_spec_iwsq(lambda, tau, wev, options)
%  Minimize squared spectral weighted barycenter weights (spectral norm) of GRQ generalized eigenvalues
% min_wsq \sum_n wev_n * wsq_n^2 * lambda_n^(-2),     (wsq_n = |w_n|^2) s.t.   
% \sum_n wsq_n = 1,   wsq_n  >= 0,    Barycentric constraints
% \sum_n wsq_n * lambda_n = tau,      GRQ constraints
% Rayleigh quotient solution: Weighted eigenvectors y = V * sqrt(wsq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% lambda:   [N x 1]     Eigenvalues
% tau:      Scalar,     Target GRQ in [min(lambda), max(lambda)] range
% wev:      [N x 1]     Channel-eigenvalue weighting

% options:  quadprog optimoptions 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% w_sq:     [N x 1]     Squared weights

arguments
    lambda (:,1) double = 1;
    tau (1,1) double = 1;
    wev (:,1) double = 1;

    options = optimoptions('quadprog', 'Display', 'none', 'MaxIterations', 10000, ...
        'OptimalityTolerance', 1e-8, 'FunctionTolerance', 1e-8, 'ConstraintTolerance', 1e-8);
end

assert(tau >= min(lambda) && tau <= max(lambda), 'tau must be in [min(lambda), max(lambda)]');
assert(numel(lambda) == numel(wev), 'wev size mismatch lambda');

N = numel(lambda);

Aeq_bary = [ones(1,N)];
beq_bary = 1;

Aeq_grq = abs(real(lambda(:)))';
beq_grq = tau;

Aeq = [Aeq_bary; Aeq_grq];
beq = [beq_bary; beq_grq];

[w_sq, fval, exitflag] = quadprog(diag(wev .* (lambda.^(-2))), zeros(N,1), [], [], Aeq, beq, zeros(N,1), ones(N,1), [], options);


