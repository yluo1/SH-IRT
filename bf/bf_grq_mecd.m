function [W, fval] = bf_grq_mecd(C, A, R, tau, options)
% Generalized Rayleigh quotient based maximum efficiency constant directivity beamformer designs

% Sample efficiency:        Spatial acoustic power average / electric power
% Constant directivity:     GRQ(A,R,W) = tau

% max_{W} W'*C*W  s.t.  W'*A*W / (W'*R*W) = tau,  W'*W = 1.

% Equivalent formulation:
% max_{W} W'*C*W  s.t.  W'*D*W = 0,  W'*W = 1,    D = A - tau * R (indefinite matrix).

% Paper reference: 
% Luo, Y. (2024, August). Constant Directivity Loudspeaker Beamforming. In 2024 32nd European Signal Processing Conference (EUSIPCO) (pp. 246-250). IEEE.

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:         [N x N] Efficiency covariance matrix
%A:         [N x N] Accept covariance matrix
%R:         [N x N] Reject covariance matrix
%tau:       Scalar, target GRQ 

%options:           struct
%options.mode:      String, compute method {'PA', 'SDP'}
%                       'PA':   Projected Ascent
%                       'SDP':  Semi-definite Program (cvx required)

% options.PA_step_size:     Step-size in direction of positive gradient
% options.PA_max_iter:      Maximum number of iterations
% options.PA_max_norm:      Maximum norm term for normalization W per step
% options.PA_func_tol:      Function difference norm tolerance for convergence
% options.PA_constr_tol:    GRQ difference from tau for satisfying constraint

%options.check_constraints: Logical, if true check for constraints satisfication of solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%W:         [N x 1] Beamformer weights
%fval:      Scalar, objective W'*C*W

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare compute methods

% rng(343);
% N = 8;
% C = randn(N) + randn(N) * 1i;     C = C * C';
% A = randn(N) + randn(N) * 1i;     A = A * A';
% R = randn(N) + randn(N) * 1i;     R = R * R';
% tau = 6;
% [W_PA, fval_PA]   = bf_grq_mecd(C, A, R, tau, 'mode', 'PA', 'check_constraints', true);
% [W_SDP, fval_SDP] = bf_grq_mecd(C, A, R, tau, 'mode', 'SDP', 'check_constraints', true);
% W_PA =  W_PA / W_PA(1);  % Normalize by first channel
% W_SDP = W_SDP / W_SDP(1);
% err = norm( W_PA - W_SDP)

arguments
    C (:,:) double {coder.mustBeComplex} = complex(1);
    A (:,:) double {coder.mustBeComplex} = complex(1);
    R (:,:) double {coder.mustBeComplex} = complex(1);
    tau (1,1) double {mustBeNonnegative} = inf;

    options.mode (1,:) char {mustBeMember(options.mode, {'PA', 'SDP'})} = 'PA';

    options.PA_step_size (1,1) double {mustBeNonnegative} = 1e-1;
    options.PA_max_iter (1,1) double {mustBeInteger, mustBePositive} = 100;
    options.PA_max_norm (1,1) double {mustBePositive} = 1e10;
    options.PA_func_tol (1,1) double {mustBeNonnegative} = 1e-8;
    options.PA_constr_tol (1,1) double {mustBeNonnegative} = 1e-8;

    options.check_constraints (1,1) logical = false;
end

N = size(C, 1);
assert(N == size(A, 1) && N == size(A, 2), 'size mismatch A, C');
assert(isequal(size(A), size(R)), 'Size A, R mismatch');
assert(N == size(C, 2), 'C must be square');

% Check tau against GRQ bounds
[W_min_max, tau_min_max] = bf_grq_ex(A, R);
if tau < tau_min_max(1) % Solution is GRQ minimizer

    W = W_min_max(:, 1);

elseif tau > tau_min_max(2) % Solution is GRQ maximizer

    W = W_min_max(:, 2);

else

    if strcmp(options.mode, 'PA') % Projected ascent
    
        W = ones(N, 1);
        fval_prev = W'*C*W;

        for i = 1:options.PA_max_iter
    
            W_step = W + options.PA_step_size * C * W;  % Move in positive objective's gradient
            W_proj = bf_grq_proj(W_step, A, R, tau);    % Project onto constraints
            W = W_proj / min( norm(W_proj), options.PA_max_norm );  % Normalize
            
            fval = W'*C*W;
            
            GRQ_curr  = real( W'*A*W / (W'*R*W) );            
            if abs(fval - fval_prev) <= options.PA_func_tol && ... % Function evaluation converged
               abs(GRQ_curr - tau) <= options.PA_constr_tol % Constraints satisfied
                break;
            else
                fval_prev = fval;
            end
            
        end
    
    elseif strcmp(options.mode, 'SDP') % Semi-definite relaxation is tight
        
        D = A - tau * R;

        cvx_clear
        cvx_begin
        variable W_mat(N,N) hermitian semidefinite
            maximize trace(W_mat * C)
            subject to
                trace(W_mat * D) == 0;
                trace(W_mat) == 1;
        cvx_end
        [W, ~] = eigs(double(W_mat), 1); % W_mat is rank-1
    
    else
        error('Unknown mode');
    end
end

% Normalize
W = W / sqrt(W'*W);
fval = real(W'*C*W);

if options.check_constraints
    tol = 1e-6;
    pass_grq_tau_c = abs(W'*A*W / (W'*R*W) - tau) <= tol
end