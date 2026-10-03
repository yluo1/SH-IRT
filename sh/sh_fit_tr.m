function [C, err, h_figs] = sh_fit_tr(X, theta, phi, max_odr, is_real, lambda, mode, options)
%Tikhonov regularization least squares: min_C ||Y(theta, phi) * C - X||^2 + C'*Q*C
%C = (Y'Y + Q)^(-1) * Y' * X

%mode = 'identity'  unity Q, penalize squared Euclidean norm of C
%Q = lambda * I

%mode = 'quad'   diagonal regularization with quadratic weighted penalty of cofficients in C w.r.t. SH degree
%Q = lambda * diag(q),   q = [0, ..., l^2, ...] for SH degree l

%mode = 'quadlin'   diagonal regularization with quadratic + linear weighted penalty of cofficients in C w.r.t. SH degree
%Q = lambda * diag(q),   q = [0, ..., l*(l+1), ...] for SH degree l

%mode = 'picard'    minimum regularization of ascending singular values below the Picard crossover index given by
%first ascending singular value normalized CDF / observation projected on left singular value normalized CDF >= lambda
%Q = V * diag(d) * V',   d = [d_*; zeros(N_C - k_*, 1)], where d_* is a minimum norm non-negative solution with linear constraints
%to ensure ascending Picard ratios and regularized singular value sequences.

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X:             [N x M] N measurements of M functions

%theta:         [N x 1] Co-latitude [0, pi]
%phi:           [N x 1] Azimuth [0, 2 * pi)

%max_odr:       Max SH order 
%is_real:       Logical, if true, evaluate real SH

%lambda:        Scalar regularization term
%mode:          String, regularization method {'identity', 'quad', 'quadlin', 'picard'}

%options:               struct
%options.W:             [], [N x 1], or [N x N] weighting matrix for weighted least squares
%                       [] uniform weighting
%                       [N x 1] Positive weight vector
%                           min_C ||diag(W.^(1/2)) * (Y(theta, phi) * C - X)||^2 + C'*Q*C
%                       [N x N] Positive definite weight matrix
%                           min_C ||chol(W) * (Y(theta, phi) * C - X)||^2 + C'*Q*C

%options.picard_log_transform:  Logical, if true, solve log-transform of linear program (numerically stable for small singular values),
%                               minimizes sum_i log_gamma_i where sigma_i^2 gamma_i = sigma_i^2 + d_i

%options.enable_disp:       Logical, if true, plot fit and Picard plot
%options.disp_picard_idx:   Function index to display for Picard plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:                 [(max_odr + 1)^2 x M] SH coefficients
%err:               Scalar, norm( W^(1/2) * (Y * C - X) )
%h_figs:            Handle to figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare Tikhonov regularization modes for random field + noise

% rng(441);
% max_odr = 6;
% is_real = false;
% N_C = (max_odr + 1)^2;
% N_pts = (max_odr + 1)^2;
% 
% %Reference field
% C_ref = sh_rand(max_odr, 1, is_real);
% 
% %Sample at spherical coordinates
% [theta, phi] = sh_grd_rand(N_pts);  %Random over sphere
% 
% X = sh_dec(C_ref, theta, phi, is_real);
% %Add noise
% X = X + (randn(size(X)) + randn(size(X)) * 1i) * 5e-2;

%lambda = 0.1;
%C_ls = sh_fit_tr(X, theta, phi, max_odr, is_real, 0);
%C_identity = sh_fit_tr(X, theta, phi, max_odr, is_real, lambda, 'identity');
%C_quad = sh_fit_tr(X, theta, phi, max_odr, is_real, lambda, 'quad');
%C_quadlin = sh_fit_tr(X, theta, phi, max_odr, is_real, lambda, 'quadlin');
%C_picard = sh_fit_tr(X, theta, phi, max_odr, is_real, 0.5, 'picard', 'enable_disp', true);

%dB_lim = [-40, 20];
%sh_plt(C_ref, 'mercator', is_real, 'title_name', 'Reference', 'dB_lim', dB_lim);
%sh_plt(C_ls, 'mercator', is_real, 'title_name', 'Non-Regularized Least Squares', 'dB_lim', dB_lim);
%sh_plt(C_identity, 'mercator', is_real, 'title_name', 'TR Least Squares identity', 'dB_lim', dB_lim);
%sh_plt(C_quad, 'mercator', is_real, 'title_name', 'TR Least Squares quaddeg', 'dB_lim', dB_lim);
%sh_plt(C_quadlin, 'mercator', is_real, 'title_name', 'TR Least Squares quadlindeg', 'dB_lim', dB_lim);


arguments
    X (:,:) double {coder.mustBeComplex} = complex(0);

    theta (:,1) double = [0];
    phi   (:,1) double = [0];

    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;

    is_real (1,1) logical = false;
    
    lambda (1,1) double {mustBeNonnegative} = 0;

    mode (1,:) char {mustBeMember(mode, {'identity', 'quad', 'quadlin', 'picard'})} = 'identity';

    options.W (:, :) double = [];
        
    options.picard_log_transform (1,1) logical = false;

    options.enable_disp (1,1) logical = false;
    options.disp_picard_idx (1,1) {mustBeNonnegative, mustBeInteger} = 1;
end

% Compute SH bases
Y = sh_val(max_odr, theta, phi, is_real);  %[N x (max_odr + 1)^2]
[N, N_C] = size(Y);
M = size(X, 2);

% Apply optional weighting matrix
if ~isempty(options.W)
    if numel(options.W) == N

        assert(all(options.W(:) > 0), 'options.W must be positive');
        W_sqrt = sqrt(options.W(:));
        Y = bsxfun(@times, Y, W_sqrt);
        X = bsxfun(@times, X, W_sqrt);

    elseif isequal(size(options.W),  [N, N])        

        [W_sqrt, flag] = chol(options.W);
        assert(flag == 0, 'options.W is not positive definite');
        Y = W_sqrt * Y;
        X = W_sqrt * X;

    else
        error('options.W size unsupported');
    end
end

% Compute regularization matrix
if any(strcmp(mode, {'identity', 'quad', 'quadlin'}))

    if strcmp(mode, 'identity')
    
        Q = lambda * eye(N_C);
        
    elseif strcmp(mode, 'quad')
    
        q = zeros(N_C, 1);
        for l = 0:max_odr
            idx = (l + 1)^2 - l;
            q(idx) = l^2;
        end
        Q = lambda * diag(q);
        
    elseif strcmp(mode, 'quadlin')
    
        q = zeros(N_C, 1);
        for l = 0:max_odr
            idx = (l + 1)^2 - l;
            q(idx) = l * (l + 1);
        end
        Q = lambda * diag(q);

    else
        error('Unsupported mode');
    end

    % Compute least squares solution
    C = (Y'*Y + Q) \ (Y' * X);

elseif strcmp(mode, 'picard')

    [U, S, V] = svd(Y, 'econ'); % Descending singular values    
    s_list = diag(S);
    s_list_ascend = flipud(s_list);
    U_ascend = fliplr(U);        
    V_ascend = fliplr(V);
     
    cdf_s = cumsum(s_list_ascend) / sum(s_list_ascend);

    C = zeros([(max_odr + 1)^2, M]);
    for m = 1:M
        U_ascend_X_m = abs(U_ascend' * X(:, m)); %[N_s_list x 1]
    
        cdf_ux = cumsum(U_ascend_X_m) / sum(U_ascend_X_m);
    
        idx_picard_cross = find((cdf_s ./ cdf_ux) >= lambda, 1, 'first');
        
        N_PC = idx_picard_cross - 1;    
        if N_PC > 0 % Setup inequality constraints
            
            A = zeros(N_PC, N_PC);
            b = zeros(N_PC, 1);

            if options.picard_log_transform % log-transform form

                log_s_list_ascend = log(s_list_ascend);
                log_U_ascend_X_m = log(U_ascend_X_m);

                for n = 1:N_PC

                    A_m = zeros(1, N_PC); % Monotonic increasing Picard ratio

                    if n < N_PC
                
                        A_m([n, n+1]) = [-1, 1];
                        b_m = log_U_ascend_X_m(n+1) - log_U_ascend_X_m(n) - log_s_list_ascend(n+1) + log_s_list_ascend(n);
    
                    else % End points
            
                        A_m(n) = [-1];
                        b_m = log_U_ascend_X_m(n+1) - log_U_ascend_X_m(n) + log_s_list_ascend(n);
                                   
                    end
        
                    A(n, :)   = A_m;
                    b(n)      = b_m;
    
                end
                   
                % Normalize
                % A_max_abs = max(abs(A), [], 2);
                % A = bsxfun(@rdivide, A, A_max_abs);
                % b = b ./ A_max_abs;
    
                % Solve
                options_linprog = optimoptions("linprog", MaxIterations=1000, ConstraintTolerance=1e-10, OptimalityTolerance=1e-10);
                [log_gamma_PC, fval, exitflag] = linprog(ones(N_PC, 1), A, b, [], [], ones(N_PC, 1), inf(N_PC, 1), options_linprog);
                
                gamma_PC = exp(log_gamma_PC);
                d_PC = s_list_ascend(1:N_PC).^2 .* (gamma_PC - 1);

            else % linear form

                for n = 1:N_PC
                    A_m = zeros(1, N_PC); % Monotonic increasing Picard ratio
                    if n < N_PC
            
                        A_m([n, n+1]) = [-U_ascend_X_m(n+1) * s_list_ascend(n+1), U_ascend_X_m(n) * s_list_ascend(n)];
                        b_m = U_ascend_X_m(n+1) * s_list_ascend(n+1) * s_list_ascend(n)^2 - U_ascend_X_m(n) * s_list_ascend(n) * s_list_ascend(n+1)^2;

                    else % End points
            
                        A_m(n) = [-U_ascend_X_m(n+1) * s_list_ascend(n+1)];
                        b_m = U_ascend_X_m(n+1) * s_list_ascend(n+1) * s_list_ascend(n)^2 - U_ascend_X_m(n) * s_list_ascend(n) * s_list_ascend(n+1)^2;                

                    end
        
                    A(n, :)   = A_m;
                    b(n)      = b_m;
                end
                % Solve
                [d_PC, fval, exitflag] = linprog(ones(N_PC, 1), A, b, [], [], zeros(N_PC, 1), inf(N_PC, 1));

            end

            d = zeros(size(V, 2), 1);
            d(1:N_PC) = d_PC;
            Q = V_ascend * diag(d) * V_ascend';
            
            if N_C > N % Augment with smallest eigenvalue paired with null-space
                V_null = null(V');
                Q = Q + V_null * min(s_list_ascend.^2 + d) * V_null'; 
            end

            % Check eigenvalues if N_C <= N
            % [V_tmp, D_tmp] = eig(Y'*Y + Q);
            % norm(sort(diag(D_tmp), 'descend') - sort(s_list_ascend.^2 + d, 'descend'))
            ;
    
        else
            Q = zeros(N_C);
        end
    
        % Compute least squares solution
        C(:, m) = (Y'*Y + Q) \ ( Y' * X(:, m) );
    end

else
    error('Unknown mode');
end

% Compute error
err = norm(Y * C - X);

h_figs = [];
if options.enable_disp && coder.target("MATLAB")
    % Plot fit
    h_figs = sh_plt(C, 'mercator', is_real, 'disp_theta_phi', [theta, phi], 'title_name', ['TR Fit ', mode, ' \lambda = ', num2str(lambda)]);

    if strcmp(mode, 'picard')

        % Picard plot 
        [U, S, V] = svd(Y, 'econ'); % Descending singular values
        s_list = diag(S);
        N_s_list = numel(s_list);

        [~, idx_descend] = sort(s_list, 'descend'); %Descending singular values
        U_descend = U(:, idx_descend);
        s_list_descend = s_list(idx_descend);
        U_descend_X = abs(U_descend' * X(:, options.disp_picard_idx)); %[N_s_list x 1]
        U_ascend_X  = flipud(U_descend_X);    

        cdf_s = cumsum(s_list_ascend) / sum(s_list_ascend);
        cdf_ux = cumsum(U_ascend_X) / sum(U_ascend_X);

        s_list_ascend_regu = (s_list_ascend.^2 + d) ./ s_list_ascend; 
        cdf_sr = cumsum(s_list_ascend_regu) / sum(s_list_ascend_regu);

        idx_picard_cross = find( cdf_s>= cdf_ux, 1, 'first');
        idx_picard_cross_trunc = find((cdf_s ./ cdf_ux) >= lambda, 1, 'first');
    
        fontsize = 18;    
        h_picard_fig = figure;
        h_picard_fig.Position = [100, 100, 800, 600];
        h_figs{end+1} = h_picard_fig;
        tiledlayout(2, 1);
        nexttile;
        semilogy(1:N_s_list, s_list_descend, 'r*-', 1:N_s_list, U_descend_X, 'bs-', ...
            1:N_s_list, U_descend_X ./ s_list_descend, 'o-', ...
            1:N_s_list, flipud(s_list_ascend_regu), 'm--', ...
            1:N_s_list, flipud(sqrt(d)), 'g-', ...
             1:N_s_list, U_descend_X ./ flipud(s_list_ascend_regu), 'k--', ...
             'linewidth', 1.5);
        grid on; axis tight;
        xlabel('Descending Singular Value Index i', 'fontsize', fontsize);
        ylabel('Magnitude', 'fontsize', fontsize);
        title('Picard Plot', 'fontsize', fontsize + 1);
        set(gca, 'fontsize', fontsize - 1);

        h_lg = legend('$\sigma_i$', '$|u_i^H x|$', '$|u_i^H x| / \sigma_i$', ...
            '$(\sigma_i^2 + d_i) / \sigma_i$', ...
            '$\sqrt{d_i}$', ...
            '$|u_i^H x| \frac{ \sigma_i} {\sigma_i^2 + d_i}$', ...
            'location' ,'best', 'interpreter', 'latex', 'NumColumns', 2); 

        set(h_lg, 'fontsize', fontsize + 1);
    
        nexttile;
        semilogy(1:N_s_list, cdf_s, 'r*-', 1:N_s_list, cdf_ux, 'bs-',  ...
            1:N_s_list, cdf_sr, 'md-', 'linewidth', 1.5); hold on;
        xline(idx_picard_cross, 'm--', 'linewidth', 2);
        xline(idx_picard_cross_trunc, 'k--', 'linewidth', 2);
       
        grid on; axis tight;
        xlabel('Ascending Singular Value Index i', 'fontsize', fontsize);
        ylabel('Magnitude', 'fontsize', fontsize);
        title('Normalized Cumulative Distribution Function', 'fontsize', fontsize + 1);
        set(gca, 'fontsize', fontsize - 1);
        if M == 1
            h_lg = legend('CDF($\sigma_i$)', 'CDF($|u_i^H x|$)',        'CDF($(\sigma_i^2 + d_i) / \sigma_i$)', 'Picard Cross: $\lambda = 1$',  ['Picard Cross: $\lambda = ', num2str(lambda), '$'], 'location' ,'best', 'interpreter', 'latex'); 
        else
            h_lg = legend('CDF($\sigma_i$)', 'CDF(mean($|u_i^H x|$))',  'CDF($(\sigma_i^2 + d_i) / \sigma_i$)', 'Picard Cross: $\lambda = 1$',  ['Picard Cross: $\lambda = ', num2str(lambda), '$'], 'location' ,'best', 'interpreter', 'latex'); 
        end
        set(h_lg, 'fontsize', fontsize - 1);
    end

end