function [W, fval, exitflag] = bf_rmvb(v, R, P, options)
% Robust minimum variance beamformer design with optional null space,
% white-noise gain, and channel gain constraints

% min_{W, lambda} W'*R*W   s.t. 
% Re(W'*v) - ||P^(1/2) * W|| >= 1,

% Reference: Lorenz, R. G., & Boyd, S. P. (2005). Robust minimum variance beamforming.
% IEEE transactions on signal processing, 53(5), 1684-1696.

% Optional constraints:
% Null space:           U'*W = 0
% White-noise gain:     W'*W <= WNG
% Channel gains:        |W(n)|^2 <= CG(n)

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%v:         [N x 1]     Steering vector for N channels
%R:         [N X N]     [N x N] Rejection covariance matrix
%P:         [N x N]     [N x N] Elipsoid covariance matrix

%options:       struct

%options.U:     Empty or [N x N_U]  Null space vectors
%               Disable constraints if empty
%options.WNG:   [], scalar, Maximum white-noise gain
%               Disable constraints if empty
%options.CG:    [], scalar, [N x 1]
%               Disable constraints if empty. Replcate across channels if scalar.

%options.check_constraints: Logical, if true check for constraints satisfication of solution

%options.options_coneprog:  Options for coneprog.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%W:         [N x 1]
%fval:      Scalar, W'*R*W
%exitflag:  1 if success, otherwise see coneprog.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Compare RMVB with MVDR

%rng(235);
%N = 3;
%v = rand(N, 1) + randn(N, 1) * 1i;
% R = randn(N) + randn(N)*1i; R = R*R';
% R = R - 0.999 * eigs(R, 1, 'smallestreal') * eye(N);
% P = randn(N) + randn(N)*1i; P = P*P'; P = 1e-3 * P;
% [W_rmvb, fval_rmvb, exitflag_rmvb] = bf_rmvb(v, R, P)

% %Compare with MVDR
%[W_mvdr, fval_mvdr] = bf_mvdr(v, R)

% %Check constraints
% W_rmvb' * v
% W_mvdr' * v

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Test RMVB optional constraints

%rng(235);
%N = 6;
%N_U = 2;
%v = rand(N, 1) + randn(N, 1) * 1i;
% R = randn(N) + randn(N)*1i; R = R*R';
% R = R - 0.999 * eigs(R, 1, 'smallestreal') * eye(N);
% P = randn(N) + randn(N)*1i; P = P*P'; P = 1e-3 * P;
% U = randn(N, N_U) + randn(N, N_U) * 1i;
% WNG = 0.4;
% CG = 0.06 * ones(N, 1);
% [W_rmvb, fval_rmvb, exitflag_rmvb] = bf_rmvb(v, R, P, 'U', U, 'WNG', WNG, 'CG', CG, 'check_constraints', true)

% test_pri_c = W_rmvb' * v
% test_u_c = U' * W_rmvb
% test_wng_c = W_rmvb' * W_rmvb
% test_cg_c = abs(W_rmvb).^2

arguments
    v (:,1) double {coder.mustBeComplex} = complex(1);
    R (:,:) double {coder.mustBeComplex} = complex(1);
    P (:,:) double {coder.mustBeComplex} = complex(1);

    options.mode_sqrt (1,:) char {mustBeMember(options.mode_sqrt, {'eig', 'chol'})} = 'eig';
    
    options.U (:,:) double {coder.mustBeComplex} = [];
    options.WNG (:,1) double = [];
    options.CG (:,1) double = [];

    options.check_constraints (1,1) logical = false;

    options.options_coneprog = optimoptions("coneprog",Display="final", ...
        MaxIterations=10000, ConstraintTolerance=1e-7, OptimalityTolerance=1e-7);  
end

N = numel(v);

assert(N == size(R, 1) && N == size(R, 2), 'v and R size mismatch');

%W_aug = [real(W); imag(W); mu]; % Variable layout

% Pre-compute matrices
if strcmp(options.mode_sqrt, 'eig')

    [V, D] = eig(R);
    R_sqrt = V * diag(sqrt(diag(D)));   R_sqrt = R_sqrt';
    [V, D] = eig(P);
    P_sqrt = V * diag(sqrt(diag(D)));   P_sqrt = P_sqrt';

elseif strcmp(options.mode_sqrt, 'chol')

    R_sqrt = chol(R);
    P_sqrt = chol(P);

else
    error('Unsupported mode_sqrt');
end

R_sqrt_re_im = [real(R_sqrt), -imag(R_sqrt); imag(R_sqrt), real(R_sqrt)];
P_sqrt_re_im = [real(P_sqrt), -imag(P_sqrt); imag(P_sqrt), real(P_sqrt)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Primary constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective as constraint
% minimize mu,  w'*R*w <= mu,   mu >= 0

Asoc = [ R_sqrt_re_im, zeros(2 * N, 1) ];
bsoc = zeros(2 * N, 1);
dsoc = [zeros(2 * N, 1); 1];
gamma = 0;
socConstraints(1) = secondordercone(Asoc,bsoc,dsoc,gamma);

% Steering vector and robustness constraints
Asoc = [ P_sqrt_re_im, zeros(2 * N, 1) ];
bsoc = zeros(2 * N, 1);
dsoc = [real(v); imag(v); 0];
gamma = 1;
socConstraints(2) = secondordercone(Asoc,bsoc,dsoc,gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optional constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Null constraints
if ~isempty(options.U)

    assert(size(options.U, 1) == N, 'options.U must be size [N x *]');
    N_U = size(options.U, 2);

    U_H = options.U'; % Conjugate transpose
    Aeq = [[real(U_H), -imag(U_H); imag(U_H), real(U_H)], zeros(2 * N_U, 1)];
    beq = zeros(2 * N_U, 1);

else
    Aeq = [];
    beq = [];
end

% White-noise gain constraints
if ~isempty(options.WNG)

    assert(isscalar(options.WNG), 'WNG must be scalar');
    Asoc = [ eye(2 * N), zeros(2 * N, 1) ];
    bsoc = zeros(2 * N, 1);
    dsoc = zeros(2 * N + 1, 1);
    gamma = -sqrt(options.WNG);
    socConstraints(end+1) = secondordercone(Asoc,bsoc,dsoc,gamma);
end

% Channel gain constraints
if ~isempty(options.CG)

    assert(isscalar(options.CG) || numel(options.CG) == N, 'CG must be scalar or size [N x 1]');

    for n = 1:N
        Asoc = zeros(2, 2 * N + 1);
        Asoc(1, n) = 1;
        Asoc(2, N + n) = 1;
        bsoc = [0; 0];
        dsoc = zeros(2 * N + 1, 1);
        if isscalar(options.CG)
            gamma = -sqrt(options.CG);
        else
            gamma = -sqrt(options.CG(n));
        end
        socConstraints(end+1) = secondordercone(Asoc,bsoc,dsoc,gamma);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = [zeros(2 * N, 1); 1];
lb = [-inf(2 * N, 1); 0]; % Lower bound

[W_aug, ~, exitflag] = coneprog(f, socConstraints, [], [], Aeq, beq, lb, [], options.options_coneprog);

if exitflag == 1

    W = W_aug(1:N) + 1i * W_aug((N+1):(2*N));
    fval = real(W'*R*W);

    % Check if solution satisfies constraints
    if options.check_constraints
    
        tol = 1e-6;
        
        pass_pri_c = real(W' * v) - sqrt(real(W'*P*W)) >= 1 - tol    
        if ~isempty(options.U)
            pass_null_c = all(abs(options.U' * W) <= tol)
        end
        if ~isempty(options.WNG)
            pass_WNG_c = real(W'*W) <= options.WNG + tol
        end
        if ~isempty(options.CG)
            pass_CG_c = all( abs(W .* conj(W)) <= options.CG + tol )
        end
    
    end

else 

    W = [];
    fval = [];

end

