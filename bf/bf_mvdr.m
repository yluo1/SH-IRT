function [W, fval, SINR, WNG] = bf_mvdr(v, R, kappa, mode)
% Minimum variance distortionless response beamformer design with condition number constraints
% min_{W, sigma} W'*R(kappa)*W   s.t.    W'*v = 1,   

% mode = 'diag'   where R is diagonal loaded R(kappa) = R + sigma * I  s.t. cond(R + sigma * I) <= kappa 
% =>    (lambda_max + sigma) / (lambda_min + sigma) <= kappa
% =>    sigma >= (lambda_max - lambda_min) / (kappa - 1)
% =>    sigma = (lambda_max - lambda_min) / (kappa - 1) at minimizer

% mode = 'thres' where  R(kappa) = V*D*V' for eigenvector matrix V s.t. cond(D) <= kappa 
% Eigenvalues below threshold are limited:
% D = diag(max(lambda, lambda_max/kappa))

% mode = 'exp'  where  R(kappa) = V*D*V' for eigenvector matrix V s.t. cond(D) <= kappa 
% Eigenvalues are exponentiated s.t. cond(D) <= kappa 
% D = diag(lambda.^p)
% Solves for max_p  s.t.  (lambda_max / lambda_min)^p <= kappa,      0 <= p <= 1
% =>    p = min(1, log(kappa) / (log(lambda_max) - log(lambda_min)) )

% mode = 'sse'  (sum-of-squared spectral error) where  R(kappa) = V*D*V' for eigenvector matrix V s.t. cond(D) <= kappa 
% D = diag(lambda_kappa);
% Minimize squared difference of eigenvalues lambda_kappa from original eigenvalues lambda s.t. 
% monotonicity and condition number contraint via quadratic programming

% mode = 'ase'  (sum-of-absolute spectral error) where  R(kappa) = V*D*V' for eigenvector matrix V s.t. cond(D) <= kappa 
% D = diag(lambda_kappa);
% Minimize absolute difference of eigenvalues lambda_kappa from original eigenvalues lambda s.t. 
% monotonicity and condition number contraint via linear programming with slack variables

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%v:         [N x 1]     Steering vector for N channels
%R:         [N X N]     Rejection covariance matrix
%kappa:                 Scalar, condition number constraint >= 1, inf = unconstrained
%                       Choose kappa <= SNR
%mode:                  String, condition number constraint method {'diag', 'thres', 'trunc', 'exp', 'sse', 'ase'}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%W:         [N x 1]     Beamformer weights (conjugated)
%fval:                  Scalar, objective:          W'*R(kappa)*W
%SINR:                  Scalar, signal-to-interference-plus-noise ratio:  |W'*v|^2 / real(W'*R*W)
%WNG:                   Scalar, white-noise gain:   W'*W

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: See tst_bf_mvdr.m

arguments
    v (:,1) double {coder.mustBeComplex} = complex(1);
    R (:,:) double {coder.mustBeComplex} = complex(1);
    kappa (1,1) double {mustBePositive} = inf;
    mode (1,:) char {mustBeMember(mode, {'diag', 'thres', 'exp', 'sse', 'ase'} )} = 'diag';
end

N = numel(v);

assert(N == size(R, 1) && N == size(R, 2), 'v and R size mismatch');
assert(kappa >= 1, 'Require condition number kappa >= 1');

if strcmp(mode, 'diag') % Shape preserving across entire eigenvalues by upwards shifting all eigenvalues
    % minimizes total difference of normalized cumulative eigenvalues and constrained eigenvalues

    lambda_max = eigs(R, 1, 'largestreal');
    lambda_min = eigs(R, 1, 'smallestreal');

    sigma = (lambda_max - lambda_min) / (kappa - 1);

    R_kappa = R + sigma * eye(N);
    R_kappa_inv_v = R_kappa \ v; 

elseif strcmp(mode, 'thres') % Shape preserving upto threshold, then brick-wall limit small eigenvalues
    % minimizes total difference between original and constrained eigenvalues

    [V, D] = eig(R);
    d = diag(D);
    [lambda, idx] = sort(d, 'descend');
    V = V(:, idx);
    
    lambda_max = lambda(1);
    
    lambda_kappa = max(lambda, lambda_max/kappa);
    D_kappa = diag(lambda_kappa);
    D_kappa_inv = diag( 1./ lambda_kappa );
    R_kappa = V * D_kappa * V';
    R_kappa_inv_v = (V * D_kappa_inv * V') * v;


elseif strcmp(mode, 'exp') % Shape smoothing eigenvalues via soft-compression

    [V, D] = eig(R);
    d = diag(D);
    [lambda, idx] = sort(d, 'descend');
    V = V(:, idx);
    
    lambda_max = lambda(1);
    lambda_min = lambda(end);
    p = min(1, log(kappa) / (log(lambda_max) - log(lambda_min)));

    lambda_kappa = lambda.^p;

    D_kappa = diag(lambda_kappa);
    D_kappa_inv = diag( lambda_kappa.^(-1) );
    R_kappa = V * D_kappa * V';
    R_kappa_inv_v = (V * D_kappa_inv * V') * v;

elseif strcmp(mode, 'sse') % Minimize squared spectral error in eigenvalues

    [V, D] = eig(R);
    d = diag(D);
    [lambda, idx] = sort(d, 'descend');
    V = V(:, idx);

    H = eye(N);
    f = -real(lambda(:));
    A = [];
    b = [];
    % Monotonicity constraints: lambda_kappa(n) >= lambda_kappa(n+1)
    for n = 2:N
        A_n = zeros(1,N); A_n(n) = 1; A_n(n-1) = -1; 
        A = [A; A_n];
        b = [b; 0];
    end
    % Condition number constraint: lambda_kappa(1) / lambda_kappa(N) <= kappa
    A_cond = zeros(1,N); A_cond(1) = 1; A_cond(N) = -kappa;
    b_cond = 0;
    A = [A; A_cond];
    b = [b; b_cond];
    
    % Non-negative eigenvalues: lambda_kappa(n) >= 0
    lb = zeros(N, 1);
    
    %Options
    options_quadprog = optimoptions('quadprog', 'Display', 'none', 'MaxIterations', 10000, ...
        'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-8, 'ConstraintTolerance', 1e-8);
    %Solve
    [lambda_kappa, fval_lambda, exitflag_lambda] = quadprog(H, f, A, b, [], [], lb, [], [], options_quadprog);
    
    D_kappa = diag(lambda_kappa);
    D_kappa_inv = diag( lambda_kappa.^(-1) );
    R_kappa = V * D_kappa * V';
    R_kappa_inv_v = (V * D_kappa_inv * V') * v;

elseif strcmp(mode, 'ase') % Minimize absolute spectral error in eigenvalues

    [V, D] = eig(R);
    d = diag(D);
    [lambda, idx] = sort(d, 'descend');
    V = V(:, idx);

    % variables order [lambda_kappa; slack] in [2 * N x 1]
    f = [zeros(N, 1); ones(N, 1)];
    A = [];
    b = [];
    % Monotonicity constraints: lambda_kappa(n) >= lambda_kappa(n+1)
    for n = 2:N
        A_n = zeros(1,2*N); A_n(n) = 1; A_n(n-1) = -1; 
        A = [A; A_n];
        b = [b; 0];
    end
    % Condition number constraint: lambda_kappa(1) / lambda_kappa(N) <= kappa
    A_cond = zeros(1,2*N); A_cond(1) = 1; A_cond(N) = -kappa;
    b_cond = 0;
    A = [A; A_cond];
    b = [b; b_cond];
    % Slack variable constraints slack(n) >= |lambda_kappa(n) - lambda(n)|
    for n = 1:N
        A_n = zeros(1,2*N); A_n(n) = -1; A_n(N+n) = -1;
        A = [A; A_n]; 
        b = [b; -real(lambda(n))];

        A_n = zeros(1,2*N); A_n(n) = 1; A_n(N+n) = -1;
        A = [A; A_n]; 
        b = [b; real(lambda(n))];
    end
    % Non-negative eigenvalues and slack: lambda_kappa(n) >= 0, slack(n) >= 0
    lb = zeros(2 * N, 1);

    %Options
    options_linprog = optimoptions("linprog", 'Display', 'none', 'MaxIterations', 1000, ...
        'ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-10);   
    %Solve
    [lambda_kappa_slack, fval_lambda, exitflag_lambda] = linprog(f, A, b, [], [], lb, [], options_linprog);
    lambda_kappa = lambda_kappa_slack(1:N);

    D_kappa = diag(lambda_kappa);
    D_kappa_inv = diag( lambda_kappa.^(-1) );
    R_kappa = V * D_kappa * V';
    R_kappa_inv_v = (V * D_kappa_inv * V') * v;

else
    error('Unknown mode');
end

% Output
W = R_kappa_inv_v / (v' * R_kappa_inv_v); % Normalize for distortionless constraint
fval = real(W'*R_kappa*W);
SINR = abs(W'*v)^2 / real(W'*R*W);
WNG = W'*W;
