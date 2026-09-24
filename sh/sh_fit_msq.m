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
%                       'SOMS':     Sum-of-magnitude squared least-squares semi-definite program (cvx)
%                       'MOMS':     Mix-of-magnitude squared least-squares semi-definite program (cvx)
%                       'MP':       Mixture power least-squares via constrained optimization (fmincon)

%options:       struct

%options.C0:                [(floor(max_odr/2) + 1)^2 x K0]  Initial guesses SH coefficients for mode = 'MS', real-valued

%options.B:                 [(floor(max_odr/2) + 1)^2 x K]  Dictionary of SH coefficients for mode = 'MOMS', 'MP'
%options.B0:                [K x K0] Initial guesses dictionary weights for mode = 'MP'

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

    options.C0 (:,:) double {mustBeReal} = [];
    options.B (:,:) double = [];
    options.B0 (:,:) double {mustBeReal}  = [];

    options.is_pdf (1,1) logical = false;

    options.options_fmincon = optimoptions("fmincon", ...
            Algorithm="interior-point", EnableFeasibilityMode=false, ...
            SpecifyObjectiveGradient=true, SpecifyConstraintGradient=true, ...
            Display="off", ScaleProblem=true, ...
            checkGradients=false, FiniteDifferenceStepSize=1e-12, FiniteDifferenceType='central', ...
            FunctionTolerance=1e-8, OptimalityTolerance=1e-10, StepTolerance=1e-8, ...
            MaxIterations=1000, MaxFunctionEvaluations = 10000);

end

Y = sh_val(max_odr, theta, phi, true);  %[N x (max_odr + 1)^2], real

[N, M] = size(X);
N_D = (max_odr + 1)^2;
N_C = (floor(max_odr/2) + 1)^2;

Y_half = Y(:, 1:N_C);
output = [];

if strcmp(mode, 'MS')
  
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
            [C_k, fval_k, exitflag_k] = fmincon(@(x) sh_msq_err(x, X(:, m), Y_half), C0(:, k), ...
                [], [], [], [], [], [], nonlcon, options.options_fmincon);
            if fval_k < fval(m)
                fval(m) = fval_k;
                C(:, m) = C_k;
                exitflag(m) = exitflag_k;
            end
        end
    end
    
    D = sh_msq(C, true); %Real
    D = [D; zeros(N_D - size(D, 1), M)]; %Real

    if ~is_real
        C = sh_re2cpx(C);
    end

elseif strcmp(mode, 'SOMS')
    
    A_mat = cell(N, 1);
    for n = 1:N
        A_mat{n} = Y_half(n,:)' * Y_half(n,:);
    end

    C = zeros([N_C, N_C, M]);
    D = zeros(N_D, M);

    for m = 1:M
        cvx_clear
        cvx_begin
            variable Q(N_C, N_C) symmetric %Q ~ C*C'
            expression f_C(N,1)
            
            for n = 1:N %Sum-of-square polynomials evaluated at theta(n), phi(n)
                f_C(n) = trace(A_mat{n} * Q);
            end

            %Least-squares objective: minimize(norm(f_C - X(:, m), 2))
            minimize(sum_square(f_C - X(:, m)))
           
            %Constraints
            subject to
                %Semi-definite cone
                Q == semidefinite(N_C);
                
                if options.is_pdf %Semi-definite relaxation
                    trace(Q) == 1;
                end

        cvx_end

        [C_eig_vec, C_eig_val] = eigs(double(Q), 1, 'largestabs'); %Largest eigenvector-eigenvalue pair

        %Check tightness
        eig_ratio = max(C_eig_val) / trace(Q);
        if eig_ratio <= 0.95
            disp(['Warning: Not tight as largest eigenvalue / trace: ', num2str(max(C_eig_val) / trace(Q))] );
        end

        %Reconstruct D as sum-of-magnitude square
        [C_eig_vec, C_eig_val] = eig(double(Q)); %Weighted eigenvector-eigenvalue pair
        C_eig_val = diag(C_eig_val);
        C(:, :, m) = C_eig_vec * diag(sqrt(C_eig_val));
                
        D_m = sum(sh_msq(C(:, :, m), true), 2); %Real
        D(1:numel(D_m), m) = D_m;

        if ~is_real
            C(:, :, m) = sh_re2cpx(C(:, :, m));
        end

    end

elseif strcmp(mode, 'MOMS')
    
    B = options.B; 
    assert(size(B, 1) == N_C, 'B size mismatch')
    K = size(B, 2);
    
    A_mat = cell(N, 1);
    for n = 1:N
        A_mat{n} = B' * Y_half(n,:)' * Y_half(n,:) * B;
    end

    C = zeros([N_C, K, M]);
    D = zeros(N_D, M);
    output.W = zeros([K, K, M]);

    for m = 1:M
        cvx_clear
        cvx_begin
            variable Q(K, K) symmetric %Q ~ W*W'
            expression f_C(N,1)
            
            for n = 1:N %Sum-of-square polynomials evaluated at theta(n), phi(n)
                f_C(n) = trace(A_mat{n} * Q);
            end

            %Least-squares objective: minimize(norm(f_C - X(:, m), 2))
            minimize(sum_square(f_C - X(:, m)))
           
            %Constraints
            subject to
                %Semi-definite cone
                Q == semidefinite(K);
                
                if options.is_pdf %Semi-definite relaxation
                    trace(Q) == 1;
                end

        cvx_end

        [W_eig_vec, W_eig_val] = eig(double(Q)); %Weighted eigenvector-eigenvalue pair
        W_eig_val = diag(W_eig_val);
        W = W_eig_vec * diag( sqrt(W_eig_val) );
        C(:, :, m) = B * W;    
        D_m = sum(sh_msq(C(:, :, m), true), 2); %Real

        D(1:numel(D_m), m) = D_m;

        if ~is_real
            C(:, :, m) = sh_re2cpx(C(:, :, m));
        end
        
        output.W(:, :, m) = W;
                
    end

elseif strcmp(mode, 'MP')

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

    % if options.is_pdf
    %     nonlcon = @(x) sh_msq_pdf_nonlcon(x);
    % else
    %     nonlcon = [];
    % end
    nonlcon = [];

    Y_half_B = Y_half * B;

    C = zeros(N_C, M);
    W = zeros(K, M);
    fval = inf(1, M);
    exitflag = zeros(1, M);

    for m = 1:M
        for k = 1:K0
            [w_k, fval_k, exitflag_k] = fmincon(@(x) sh_msq_err(x, X(:, m), Y_half_B), B0(:, k), ...
                [], [], [], [], [], [], nonlcon, options.options_fmincon);
            if fval_k < fval(m)
                fval(m) = fval_k;
                W(:, m) = w_k;
                C(:, m) = B * w_k;
                exitflag(m) = exitflag_k;
            end
        end
    end
    
    D = sh_msq(C, true); %Real
    D = [D; zeros(N_D - size(D, 1), M)]; %Real

    if ~is_real
        C = sh_re2cpx(C);
    end

    output.W = W;

else
    error('Unknown mode');
end

if ~is_real
    D = sh_re2cpx(D);
end

%Compute error
err = norm(Y * D - X);
