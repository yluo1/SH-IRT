function [W, fval] = bf_grq_proj(w, A, R, tau, mode)
% Steering vector projection onto generalized Rayleigh quotient (GRQ) equality constraints:
% Find nearest beamformer weights W = w + v from steering vector w that satisfies GRQ(A,R) = tau

% min_{v} v'*v s.t. W'*A*W / (W'*R*W) = tau,    W = w + v

% Equivalent formulation:
% min_{v} w'*w s.t. W'*D*W = 0, D = A - tau * R (indefinite matrix),    W = w + v

% Notes: tau is thresholded to be within tau_min_max (see bf_grq_ex.m)

% Paper reference: 
% Luo, Y. (2024, August). Constant Directivity Loudspeaker Beamforming. In 2024 32nd European Signal Processing Conference (EUSIPCO) (pp. 246-250). IEEE.

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%w:         [N x 1] Steering vector for N channels   
%A:         [N x N] Accept covariance matrix
%R:         [N x N] Reject covariance matrix
%tau:       Scalar, target GRQ 
%mode:      String, compute method {'secular', 'fmincon'}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%W:         [N x 1] Nearest beamformer weights W = w + v satisfying GRQ constraints
%                   scaled to have unit norm
%fval:      Scalar, squared norm objective v'*v

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compute nearest constant GRQ beamformer coefficients

% rng(343);
% N = 8;
% w = randn(N, 1) + randn(N, 1) * 1i;
% A = randn(N) + randn(N) * 1i;     A = A * A';
% R = randn(N) + randn(N) * 1i;     R = R * R';
% tau = 6;
% [W_sec, fval_sec] = bf_grq_proj(w, A, R, tau, 'secular')
% [W_fmc, fval_fmc] = bf_grq_proj(w, A, R, tau, 'fmincon')

% %Check constraints
% abs(W_sec'*A*W_sec / (W_sec'*R*W_sec) - tau) <= 1e-8
% abs(W_fmc'*A*W_fmc / (W_fmc'*R*W_fmc) - tau) <= 1e-8

arguments
    w (:,1) double {coder.mustBeComplex} = complex(1);
    A (:,:) double {coder.mustBeComplex} = complex(1);
    R (:,:) double {coder.mustBeComplex} = complex(1);
    tau (1,1) double {mustBeNonnegative} = inf;
 
    mode (1,:) char {mustBeMember(mode, {'secular', 'fmincon'})} = 'secular';
end

N = numel(w);
assert(N == size(A, 1) && N == size(A, 2), 'size mismatch v, A');
assert(isequal(size(A), size(R)), 'Size A, R mismatch');

% Check tau against GRQ bounds
[W_min_max, tau_min_max] = bf_grq_ex(A, R);
if tau < tau_min_max(1) % Solution is GRQ minimizer

    W = W_min_max(:, 1);

elseif tau > tau_min_max(2) % Solution is GRQ maximizer

    W = W_min_max(:, 2);

else

    if strcmp(mode, 'secular')

        D = A - tau * R;
        [V, E] = eig(D);

        u = V' * w;
        a = u .* conj(u);
        b = 1 ./ real(diag(E));

        lambda = sc_quad_proj(a, b, 'NewtonRaphson');
        lambda_re = lambda(imag(lambda) == 0);        
        lambda_crit = lambda_re(max(b(b < 0)) < lambda_re & lambda_re < min(b(b > 0)) );
        assert(~isempty(lambda_crit), 'Critical point not found');
        v = lambda_crit(1) * ((eye(N) - lambda_crit(1) * D) \ (D * w));

        fval = v'*v;
        W = w + v;

    elseif strcmp(mode, 'fmincon')

        options_fmincon = optimoptions('fmincon', 'FunctionTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'TolX', 1e-10, 'ConstraintTolerance', 1e-10, ...
                'MaxIter', 1000,  'Display', 'none', ...
                'SpecifyObjectiveGradient',true, 'SpecifyConstraintGradient', true,'CheckGradients', false);
    
        func = @(x) func_obj(x);

        if isreal(A) && isreal(R) && isreal(w)
        
            D = A - tau * R;
            nonlcon = @(x) f_nonlcon(x, D, w);
            x0 = ones(N, 1);
            [v, ~, exitflag]  = fmincon(func, x0, [],[],[],[],[],[], nonlcon, options_fmincon);
            assert(exitflag > 0, 'fmincon not converge');

        else

            A_re_im = [real(A) -imag(A); imag(A) real(A)];
            R_re_im = [real(R) -imag(R); imag(R) real(R)];
            w_re_im = [real(w); imag(w)];
            
            D_re_im = A_re_im - tau * R_re_im;
            nonlcon = @(x) f_nonlcon(x, D_re_im, w_re_im);
            x0 = ones(2 * N, 1);
            [v_re_im, ~, exitflag]  = fmincon(func, x0, [],[],[],[],[],[], nonlcon, options_fmincon);
            assert(exitflag > 0, 'fmincon not converge');

            v = v_re_im(1:N) + 1i * v_re_im((N+1):end);

        end
       
        fval = v'*v;
        W = w + v;

    else
        error('Unknown mode');
    end

    % Normalize
    W = W / sqrt(W'*W);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function
function [f, g] = func_obj(v)
f = v'*v;
if nargout > 1
    g = 2 * v;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-linear constraints
% Output
% cons:          [N_cons x 1]
% conseq:        [N_conseq x 1]
% gradc:         [N_vars x N_cons]
% gradceq:       [N_vars x N_conseq]
function [cons, conseq, gradc, gradceq] = f_nonlcon(v, D, w)

W = w + v;
cons = [];
conseq = W'* D * W;

if nargout > 2
   gradc = [];
   gradceq = 2 * D * W;
end