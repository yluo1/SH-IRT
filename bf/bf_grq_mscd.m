function [W, fval] = bf_grq_mscd(c, A, R, tau, options)
% Generalized Rayleigh quotient based maximum sensitivity constant directivity beamformer designs:

% Sample sensitivity:       Spatial acoustic on-axis power / electric power
% Constant directivity:     GRQ(A,R,W) = tau

% min_{W} W'*W  s.t.  W'*A*W / (W'*R*W) = tau,  W'*c = 1.

% Equivalent formulation:
% min_{W} W'*W  s.t.  W'*D*W = 0,  W'*c = 1,    D = A - tau * R (indefinite matrix).

% Paper reference: 
% Luo, Y. (2024, August). Constant Directivity Loudspeaker Beamforming. In 2024 32nd European Signal Processing Conference (EUSIPCO) (pp. 246-250). IEEE.

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%c:         [N x 1] Steering vector for N channels
%A:         [N x N] Accept covariance matrix
%R:         [N x N] Reject covariance matrix
%tau:       Scalar, target GRQ 

%options:           struct
%options.mode:      String, compute method {'secular', 'fmincon'}
%                       'secular':   Secular quadratic equation critical point
%                       'fmincon':   Constrained minimization

%options.check_constraints: Logical, if true check for constraints satisfication of solution
%options.options_fmincon:   options struct for options.mode = 'fmincon'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%W:         [N x 1] Beamformer weights
%fval:      Scalar, objective W'*W

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare compute methods

% rng(343);
% N = 8;
% c = randn(N,1) + randn(N,1) * 1i;
% A = randn(N) + randn(N) * 1i;     A = A * A';
% R = randn(N) + randn(N) * 1i;     R = R * R';
% tau = 6;
% [W_sec, fval_sec] = bf_grq_mscd(c, A, R, tau, 'mode', 'secular', 'check_constraints', true);
% [W_fmc, fval_fmc] = bf_grq_mscd(c, A, R, tau, 'mode', 'fmincon', 'check_constraints', true);
% err = norm( W_sec - W_fmc)

arguments
    c (:,1) double {coder.mustBeComplex} = complex(1);
    A (:,:) double {coder.mustBeComplex} = complex(1);
    R (:,:) double {coder.mustBeComplex} = complex(1);
    tau (1,1) double {mustBeNonnegative} = inf;

    options.mode (1,:) char {mustBeMember(options.mode, {'secular', 'fmincon'})} = 'secular';
    options.check_constraints (1,1) logical = false;

    options.options_fmincon = optimoptions('fmincon', 'FunctionTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'TolX', 1e-10, 'ConstraintTolerance', 1e-10, ...
            'MaxIter', 1000,  'Display', 'none', ...
            'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'CheckGradients', false);   
end

N = numel(c);
assert(N == size(A, 1) && N == size(A, 2), 'c and A size mismatch');
assert(N == size(R, 1) && N == size(R, 2), 'c and R size mismatch');

% Check tau against GRQ bounds
[W_min_max, tau_min_max] = bf_grq_ex(A, R);
if tau < tau_min_max(1) % Solution is GRQ minimizer

    W = W_min_max(:, 1);

elseif tau > tau_min_max(2) % Solution is GRQ maximizer

    W = W_min_max(:, 2);

else

    if strcmp(options.mode, 'secular')
    
        D = A - tau * R;
        [V, E] = eig(D);

        u = V' * c;
        a = u .* conj(u);
        b = 1 ./ real(diag(E));

        lambda = sc_quad_proj(a, b, 'NewtonRaphson');
        lambda_re = lambda(imag(lambda) == 0);        
        lambda_crit = lambda_re(max(b(b < 0)) < lambda_re & lambda_re < min(b(b > 0)) );
        assert(~isempty(lambda_crit), 'Critical point not found');
        W = lambda_crit(1) * ((eye(N) - lambda_crit(1) * D) \ c);
        W = W / (W'*c);

    elseif strcmp(options.mode, 'fmincon')

        func = @(x) func_obj(x);

        if isreal(A) && isreal(R) && isreal(w)
        
            D = A - tau * R;
            nonlcon = @(x) f_nonlcon(x, D, c);
            x0 = ones(N, 1);
            [W, ~, exitflag]  = fmincon(func, x0, [],[],[],[],[],[], nonlcon, options.options_fmincon);
            assert(exitflag > 0, 'fmincon not converge');

        else

            A_re_im = [real(A) -imag(A); imag(A) real(A)];
            R_re_im = [real(R) -imag(R); imag(R) real(R)];
            c_re_im = [real(c); imag(c)];
            
            D_re_im = A_re_im - tau * R_re_im;
            nonlcon = @(x) f_nonlcon(x, D_re_im, c_re_im);
            x0 = ones(2 * N, 1);
            [W_re_im, ~, exitflag]  = fmincon(func, x0, [],[],[],[],[],[], nonlcon, options.options_fmincon);
            assert(exitflag > 0, 'fmincon not converge');

            W = W_re_im(1:N) + 1i * W_re_im((N+1):end);

        end

    else
        error('Unknown mode');
    end

end

fval = real(W'*W);

if options.check_constraints
    tol = 1e-6;
    pass_grq_tau_c = abs(W'*A*W / (W'*R*W) - tau) <= tol
    pass_lin_eq_c = abs(W'*c - 1) <= tol
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function
function [f, g] = func_obj(W)
f = W'*W;

if nargout > 1
    g = 2 * W;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-linear constraints
% Output
% cons:          [N_cons x 1]
% conseq:        [N_conseq x 1]
% gradc:         [N_vars x N_cons]
% gradceq:       [N_vars x N_conseq]
function [cons, conseq, gradc, gradceq] = f_nonlcon(W, D, c)

cons = [];
conseq = [W'* D * W, W'*c - 1];

if nargout > 2
   gradc = [];
   gradceq = [2 * D * W, c];
end
