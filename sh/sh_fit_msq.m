function [D, C, err, output] = sh_fit_msq(X, theta, phi, max_odr, is_real, mode, options)
%Least-squares fit of magnitude squared SH expansion variants to real observations:

%\min_D \sum_{\theta, phi} |Y(theta, phi) * D - X(theta, phi)|^2    s.t.

%Magnitude square form for unknown C (mode = 'MS')
%Y(theta, phi) * D = abs(Y(theta, phi) * C)^2

%Sum-of-magnitude square form for unknown C_n (mode = 'SOMS')
%Y(theta, phi) * D = \sum_n abs(Y(theta, phi) * C_n)^2 

%Mix-of-magnitude square form for unknown w_n, and dictionary of coefficients B (mode = 'MOMS')
%Y(theta, phi) * D = \sum_n abs(Y(theta, phi) * B * w_n)^2 

%Mixture power form for unknown w, and dictionary of coefficients B (mode = 'MP')

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X:             [N x M] N measurements of M functions, real

%theta:         [N x 1] Co-latitude [0, pi]
%phi:           [N x 1] Azimuth [0, 2 * pi)

%max_odr:       Max SH order 

%is_real:       Logical, if true, evaluate real SH

%mode:          String, fitting method {'MS', 'SOMS', 'MOMS', 'MP'}
%                       'MS':       Magnitude squared least-squares via constrained optimization (fmincon)
%                       'MP':       Mixture power least-squares via constrained optimization (fmincon)
%                       'SOMS':     Sum-of-magnitude squared least-squares semi-definite program (cvx)
%                       'MOMS':     Mix-of-magnitude squared least-squares semi-definite program (cvx)

%options:       struct

%options.C0:                [(floor(max_odr/2) + 1)^2 x K0]  Initial guesses SH coefficients for mode = 'MS'
%options.B0:                [K x K0] Initial guesses dictionary weights for mode = 'MP'

%options.B:                 [(floor(max_odr/2) + 1)^2 x K]  Dictionary of SH coefficients for mode = 'MOMS', 'MP'

%options.is_pdf:            Logical, if true, enforce unity total energy constraint on D

%options.options_fmincon:   Options struct for fmincon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%D:                 [(max_odr + 1)^2 x M] SH coefficients

%C:                 [(floor(max_odr/2) + 1)^2 x M] SH coefficients for mode = 'MS'
%                   [(floor(max_odr/2) + 1)^2 x (floor(max_odr/2) + 1)^2 x M] SH coefficients for mode = 'SOMS'
%                   [(floor(max_odr/2) + 1)^2 x K x M] SH coefficients for mode = 'MOMS'

%err:               Sum-of-squared errors:   sum((Y(theta, phi) * D - X(theta, phi))^2)

%output:            struct, additional outputs
%output.W:          [K x K x M] mixture weights for mode = 'MOMS'
%                   [K x M] mixture weights for mode = 'MP'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: See tst_sh_fit_msq.m

arguments
    X (:,:) double {mustBeReal} = 1;

    theta (:,1) double = [0];
    phi   (:,1) double = [0];

    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;

    is_real (1,1) logical = false;

    mode (1,:) char {mustBeMember(mode, {'MS', 'SOMS', 'MOMS', 'MP'})} = 'MS';

    options.C0 (:,:) double = [];
    options.B (:,:) double = [];
    options.B0 (:,:) double  = [];

    options.is_pdf (1,1) logical = false;

    options.options_fmincon = optimoptions("fmincon", ...
            Algorithm="interior-point", EnableFeasibilityMode=false, ...
            SpecifyObjectiveGradient=true, SpecifyConstraintGradient=true, ...
            checkGradients=false, FiniteDifferenceStepSize=1e-12, FiniteDifferenceType='central', ...
            Display="off", ScaleProblem=true, ...
            FunctionTolerance=1e-8, OptimalityTolerance=1e-10, StepTolerance=1e-8, ...
            MaxIterations=1000, MaxFunctionEvaluations = 10000);

end

[N, M] = size(X);
N_D = (max_odr + 1)^2;
N_C = (floor(max_odr/2) + 1)^2;

output = [];

if strcmp(mode, 'MS')
  
    Y = sh_val(max_odr, theta, phi, is_real);  %[N x (max_odr + 1)^2]
    Y_half = Y(:, 1:N_C);

    if isempty(options.C0) %Initial guess   
        C0 = sh_enc_uni(floor(max_odr/2)); 
    else
        C0 = options.C0;   
    end
    assert(size(C0, 1) == N_C, 'C0 size mismatch')
    K0 = size(C0, 2);

    if options.is_pdf
        nonlcon = @(x) sh_msq_pdf_nonlcon(x);
    else
        nonlcon = [];
    end

    C = zeros(N_C, M);    
    fval = inf(1, M);
    exitflag = zeros(1, M);

    for m = 1:M
        for k = 1:K0

            if is_real %Real case

                [C_k, fval_k, exitflag_k] = fmincon(@(x) sh_msq_err(x, X(:, m), Y_half), C0(:, k), ...
                    [], [], [], [], [], [], nonlcon, options.options_fmincon);

            else % Complex case, interleave real and imag

                C0_re_im_k = reshape( [real(C0(:, k)), imag(C0(:, k))]', [], 1 );            
                Y_half_re_im = reshape( permute(cat(3, real(Y_half), imag(Y_half)), [1 3 2] ), [N, 2 * N_C]);

                [C_re_im_k, fval_k, exitflag_k] = fmincon(@(x_re_im) sh_msq_err_complex(x_re_im, X(:, m), Y_half_re_im), C0_re_im_k, ...
                    [], [], [], [], [], [], nonlcon, options.options_fmincon);

                C_k = C_re_im_k(1:2:end) + 1i * C_re_im_k(2:2:end);

            end
            if fval_k < fval(m)
                fval(m) = fval_k;
                C(:, m) = C_k;
                exitflag(m) = exitflag_k;
            end
        end
    end
    
    D_msq = sh_msq(C, is_real);
    D = [D_msq; zeros(N_D - size(D_msq, 1), M)];

elseif strcmp(mode, 'MP')

    Y = sh_val(max_odr, theta, phi, is_real);  %[N x (max_odr + 1)^2]
    Y_half = Y(:, 1:N_C);

    B = options.B; 
    assert(size(B, 1) == N_C, 'B size mismatch')
    K = size(B, 2);

    if isempty(options.B0) %Initial guess   
        B0 = zeros(K, 1); 
    else
        B0 = options.B0;   
    end
    assert(size(B0, 1) == K, 'B0 size mismatch')
    K0 = size(B0, 2);

    if options.is_pdf
        if is_real
            nonlcon = @(x) sh_msq_real_mix_pdf_nonlcon(x, B);
        else
            B_re_im = reshape( permute(cat(3, real(B), imag(B)), [3 1 2] ), [2 * N_C, K]);
            nonlcon = @(x) sh_msq_complex_mix_pdf_nonlcon(x, B_re_im);
        end
    else
        nonlcon = [];
    end

    Y_half_B = Y_half * B; % [N x K]

    C = zeros(N_C, M);
    W = zeros(K, M);
    fval = inf(1, M);
    exitflag = zeros(1, M);

    for m = 1:M
        for k = 1:K0
            if is_real %Real case

                [w_k, fval_k, exitflag_k] = fmincon(@(x) sh_msq_err(x, X(:, m), Y_half_B), B0(:, k), ...
                    [], [], [], [], [], [], nonlcon, options.options_fmincon);

            else % Complex case, interleave real and imag

                B0_re_im_k = reshape( [real(B0(:, k)), imag(B0(:, k))]', [], 1 );    
                Y_half_B_re_im = reshape( permute(cat(3, real(Y_half_B), imag(Y_half_B)), [1 3 2] ), [N, 2 * K]);
                
                [w_re_im_k, fval_k, exitflag_k] = fmincon(@(x) sh_msq_err_complex(x, X(:, m), Y_half_B_re_im), B0_re_im_k, ...
                    [], [], [], [], [], [], nonlcon, options.options_fmincon);

                w_k = w_re_im_k(1:2:end) + 1i * w_re_im_k(2:2:end);
                
            end
            
            if fval_k < fval(m)
                fval(m) = fval_k;
                W(:, m) = w_k;
                C(:, m) = B * w_k;
                exitflag(m) = exitflag_k;
            end
        end
    end
    
    D = sh_msq(C, is_real);
    D = [D; zeros(N_D - size(D, 1), M)];

    output.W = W;


elseif strcmp(mode, 'SOMS')
    
    Y = sh_val(max_odr, theta, phi, is_real);  %[N x (max_odr + 1)^2]
    Y_half = Y(:, 1:N_C);

    A_mat = cell(N, 1);
    for n = 1:N
        A_mat{n} = Y_half(n,:)' * Y_half(n,:);
        A_mat{n} = (A_mat{n} + A_mat{n}') / 2;
    end

    C = zeros([N_C, N_C, M]);
    D = zeros(N_D, M);

    for m = 1:M
        cvx_clear
        cvx_begin
            if is_real
                variable Q(N_C, N_C) symmetric semidefinite %Q ~ C*C'
            else
                variable Q(N_C, N_C) hermitian semidefinite %Q ~ C*C'
            end
            expression f_C(N,1)
            
            for n = 1:N %Sum-of-square polynomials evaluated at theta(n), phi(n)
                f_C(n) = trace(A_mat{n} * Q);
            end

            %Least-squares objective: minimize(norm(f_C - X(:, m), 2))
            minimize(sum_square(f_C - X(:, m)))
            %minimize(norm(f_C - X(:, m), 2))
           
            %Constraints
            subject to
                
                if options.is_pdf %Semi-definite relaxation
                    trace(Q) == 1;
                end

        cvx_end

        %Check tightness
        % [C_eig_vec, C_eig_val] = eigs(double(Q), 1, 'largestabs'); %Largest eigenvector-eigenvalue pair
        % eig_ratio = max(C_eig_val) / trace(Q);
        % if eig_ratio <= 0.95
        %     disp(['Warning: Not tight as largest eigenvalue / trace: ', num2str(max(C_eig_val) / trace(Q))] );
        % end

        %Reconstruct D as sum-of-magnitude square
        [C_eig_vec, C_eig_val] = eig(double(Q)); %Weighted eigenvector-eigenvalue pair
        C_eig_val = diag(C_eig_val);
        [~, idx_sort] = sort(C_eig_val, 'descend');
        C(:, :, m) = C_eig_vec(:, idx_sort) * diag(sqrt(C_eig_val(idx_sort)));
        D_m = sum(sh_msq(C(:, :, m), is_real), 2);
        D(1:numel(D_m), m) = D_m;

    end

elseif strcmp(mode, 'MOMS')
    
    Y = sh_val(max_odr, theta, phi, is_real);  %[N x (max_odr + 1)^2]
    Y_half = Y(:, 1:N_C);

    B = options.B; 
    assert(size(B, 1) == N_C, 'B size mismatch')
    K = size(B, 2);
    
    A_mat = cell(N, 1);
    for n = 1:N
        A_mat{n} = B' * Y_half(n,:)' * Y_half(n,:) * B;
        A_mat{n} = (A_mat{n} + A_mat{n}') / 2;
    end

    C = zeros([N_C, K, M]);
    D = zeros(N_D, M);
    output.W = zeros([K, K, M]);

    for m = 1:M
        cvx_clear
        cvx_begin
            if is_real
                variable Q(K, K) symmetric semidefinite %Q ~ W*W'
            else
                variable Q(K, K) hermitian semidefinite %Q ~ W*W'
            end
            expression f_C(N,1)
            
            for n = 1:N %Sum-of-square polynomials evaluated at theta(n), phi(n)
                f_C(n) = trace(A_mat{n} * Q);
            end

            %Least-squares objective: minimize(norm(f_C - X(:, m), 2))
            minimize(sum_square(f_C - X(:, m)))
            %minimize(norm(f_C - X(:, m), 2))
           
            %Constraints
            subject to
                
                if options.is_pdf %Semi-definite relaxation
                    trace(Q * (B' * B)) == 1;
                end

        cvx_end

        [W_eig_vec, W_eig_val] = eig(double(Q)); %Weighted eigenvector-eigenvalue pair
        W_eig_val = diag(W_eig_val);
        [~, idx_sort] = sort(W_eig_val, 'descend');
        W = W_eig_vec(:, idx_sort) * diag( sqrt(W_eig_val(idx_sort)) );
        C(:, :, m) = B * W;
        D_m = sum(sh_msq(C(:, :, m), is_real), 2);
        D(1:numel(D_m), m) = D_m;
           
        output.W(:, :, m) = W;
                
    end

else
    error('Unknown mode');
end

%Compute error
err = norm(Y * D - X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Complex variant with interleaved real and imag variables

%Input

%C_re_im: [2 * (P + 1)^2 x 1]       Spherical harmonic expansion coefficients real, imag interleaved
%X: [N x 1]                         N target observations, non-negative
%Y_re_im: [N x 2 * (P + 1)^2]       Spherical harmonic bases evaluations real, imag interleaved

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%f:             Scalar, sum of square errors
%g_re_im:       [2 * (P + 1)^2 x 1] gradient g(C) = df(C)/dC, real, imag interleaved

function [f, g_re_im] = sh_msq_err_complex(C_re_im, X, Y_re_im)

C = C_re_im(1:2:end) + 1i * C_re_im(2:2:end);
Y = Y_re_im(:, 1:2:end) + 1i * Y_re_im(:, 2:2:end);

[f, g] = sh_msq_err(C, X, Y);

g_re_im = reshape([real(g), imag(g)]', [], 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ineqnonlin, eqnonlin, gradineqnonlin, gradeqnonlin] = sh_msq_real_mix_pdf_nonlcon(w, B, relax_eq_constr)
%Equality constraints for unity integration of magnitude squared real spherical harmonic
%expansion density function mixtures

%C = B * w;

%int abs(Y(\omega) * C)^2 d\omega = 1 
%<=>
%C'*C = 1

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%w: [K x 1]                 Real mixture weights
%B: [(P + 1)^2 x K]         Real spherical harmonic expansion dictionary

%relax_eq_constr:           Logical, if true, relax equality constraint C'*C = 1 
%                           to inequality constraints C'*C <= 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%ineqnonlin:        Scalar or empty
%eqnonlin:          Scalar or empty
%gradineqnonlin:    [K x 1], or empty
%gradeqnonlin:      [K x 1], or empty

arguments
     w (:,1) double {mustBeReal} = 0;     
     B (:,:) double {mustBeReal} = 0;
     relax_eq_constr (1,1) logical = false;
end

C = B * w;

if relax_eq_constr
    ineqnonlin = C'*C - 1;
    eqnonlin =  [];

    if nargout > 2
        gradineqnonlin = 2 * B' * C;
        gradeqnonlin = [];
    end

else
    ineqnonlin = [];
    eqnonlin =  C'*C - 1;

    if nargout > 2
        gradineqnonlin = [];
        gradeqnonlin = 2 * B' * C;
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ineqnonlin, eqnonlin, gradineqnonlin_re_im, gradeqnonlin_re_im] = ...
    sh_msq_complex_mix_pdf_nonlcon(w_re_im, B_re_im, relax_eq_constr)
%Equality constraints for unity integration of magnitude squared complex spherical harmonic
%expansion density function mixtures

%C = B * w;

%int abs(Y(\omega) * C)^2 d\omega = 1 
%<=>
%C'*C = 1

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%w_re_im: [2 * K x 1]                 Real, imag interleaved mixture weights
%B_re_im: [2 * (P + 1)^2 x K]         Real, imag interleaved spherical harmonic expansion dictionary

%relax_eq_constr:           Logical, if true, relax equality constraint C'*C = 1 
%                           to inequality constraints C'*C <= 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%ineqnonlin:        Scalar or empty
%eqnonlin:          Scalar or empty
%gradineqnonlin_re_im:    [2 * K x 1], or empty
%gradeqnonlin_re_im:      [2 * K x 1], or empty

arguments
     w_re_im (:,1) double = 0;     
     B_re_im (:,:) double = 0;
     relax_eq_constr (1,1) logical = false;
end

w = w_re_im(1:2:end) + 1i * w_re_im(2:2:end);
B = B_re_im(1:2:end, :) + 1i * B_re_im(2:2:end, :);

C = B * w;

if relax_eq_constr
    ineqnonlin = real(C'*C) - 1;
    eqnonlin =  [];

    if nargout > 2
        gradineqnonlin = 2 * B' * C;
        gradineqnonlin_re_im = reshape([real(gradineqnonlin), imag(gradineqnonlin)]', [], 1);
        gradeqnonlin_re_im = [];
    end

else
    ineqnonlin = [];
    eqnonlin =  real(C'*C) - 1;

    if nargout > 2
        gradineqnonlin_re_im = [];
        gradeqnonlin = 2 * B' * C;
        gradeqnonlin_re_im = reshape([real(gradeqnonlin), imag(gradeqnonlin)]', [], 1);
    end

end


