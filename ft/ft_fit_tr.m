function [h, err, h_figs] = ft_fit_tr(X, w, num_taps, lambda, mode, options)
%Regularized least squares fit of FIR filters h to frequency response matrix X
%min_h ||F(w) * h - X||^2 + h'*Q*h,   h is real, F is Fourier transform matrix at w

%mode = 'identity'  unity Q, penalize squared Euclidean norm of h
%Q = lambda * I

%mode == 'gauss'    Complementary Gaussian shaped weights for regularization of h 
%Q = lambda * diag(q),     q(n) = 1 - exp( - ( (n-1) - gauss_mu )^2 / (2 * gauss_std^2) )

%mode == 'gauss'    Complementary circular Gaussian shaped weights for regularization of h 
%Q = lambda * diag(q),     q(n) = c - sum_k exp( - ( (n-1) - gauss_mu + k * num_taps )^2 / (2 * gauss_std^2) ), k = -inf to inf, c = max(exp_sum)

%mode = 'picard'    minimum regularization of ascending singular values below the Picard crossover index given by
%first ascending singular value normalized CDF / observation projected on left singular value normalized CDF >= lambda
%Q = V * diag(d) * V',   d = [d_*; zeros(N_C - k_*, 1)], where d_* is a minimum norm non-negative solution with linear constraints
%to ensure ascending Picard ratios and regularized singular value sequences.

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X:             [N x M] N frequency responses, M functions
%w:             [N x 1] Angular frequency (radians / sample) [0 to 2 * pi]

%num_taps:      Number of fitted FIR taps
%lambda:        Scalar regularization term
%mode:          String, regularization method {'identity', 'gauss', 'circulargauss', 'picard'}

%options:               Struct

%options.W:             [], [N x 1], or [N x N] Frequency weighting matrix for weighted least squares
%                       [] uniform weighting
%                       [N x 1] Positive weight vector
%                           min_h ||diag(W.^(1/2)) * (F(w) * h - X)||^2 + h'*Q*h
%                       [N x N] Positive definite weight matrix
%                           min_h ||chol(W) * (F(w) * h - X)||^2 + h'*Q*h

%options.gauss_mu:      Mean delay regularization weights for mode is 'gauss', 'circulargauss'
%options.gauss_std:     Standard deviation delay regularization weights for mode = 'gauss', 'circulargauss'

%options.picard_log_transform:  Logical, if true, solve log-transform of linear program (numerically stable for small singular values),
%                               minimizes sum_i log_gamma_i where sigma_i^2 gamma_i = sigma_i^2 + d_i

%options.enable_disp:       Logical, if true, plot filter and fitting error
%options.disp_picard_idx:   Function index to display for Picard plot
%options.Fs:                Sampling rate for display

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h:             [num_taps x M] FIR filter
%err:           Scalar, norm( W^(1/2) * (F * h - X) )
%h_figs:        Handle to figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Fit FIR to log-frequency sampled spherical piston model

% Fs = 48000;
% hz = logspace(log10(10), log10(Fs / 2), 512 + 1);
% w  = hz / Fs * (2 * pi);
% r = 3; % Evaluation distance
% radius = 0.25; % Spherical baffle radius
% C_ff = sh_enc_pist_sphere(40, pi/2, 0, radius, deg2rad(20), hz, r, 'enable_disp', true);
% X = sh_dec(C_ff, pi/2, deg2rad(15), false).';

% N_taps = 1024;
% [h_gauss, err_gauss] = ft_fit_tr(X, w, N_taps, 1e-3, 'gauss', 'gauss_mu', (r - radius) /343 * Fs, 'gauss_std', 200, 'enable_disp', true); err_gauss
% [h_id, err_id] = ft_fit_tr(X, w, N_taps, 1e-1, 'identity', 'enable_disp', true); err_id
% [h_picard, err_picard] = ft_fit_tr(X, w, N_taps, 0.5, 'picard', 'enable_disp', true); err_picard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Fit FIR to log-frequency sampled normalized spherical piston model

% Fs = 48000;
% hz = logspace(log10(10), log10(Fs / 2), 512 + 1);
% w  = hz / Fs * (2 * pi);
% r = 3; % Evaluation distance
% radius = 0.25; % Spherical baffle radius
% C_ff = sh_enc_pist_sphere(40, pi/2, 0, radius, deg2rad(20), hz, r, 'normalize_piston_axis', true, 'enable_disp', true);
% X = sh_dec(C_ff, pi/2, deg2rad(15), false).';

% N_taps = 64;
% [h_picard_log, err_picard_log] = ft_fit_tr(X, w, N_taps, 0.5, 'picard', 'picard_log_transform', true, 'enable_disp', true); err_picard_log
% [h_picard_lin, err_picard_lin] = ft_fit_tr(X, w, N_taps, 0.5, 'picard', 'picard_log_transform', false, 'enable_disp', true); err_picard_lin

% [h_id, err_id] = ft_fit_tr(X, w, N_taps, 1e-1, 'identity', 'enable_disp', true); err_id
% [h_gauss, err_gauss] = ft_fit_tr(X, w, N_taps, 1e-1, 'gauss', 'gauss_std', 20, 'enable_disp', true); err_gauss
% [h_cgauss, err_cgauss] = ft_fit_tr(X, w, N_taps, 1e-1, 'circulargauss', 'gauss_std', 20, 'enable_disp', true); err_cgauss

arguments
    X (:,:) double = 1;
    w (:,1) double = 0;

    num_taps (1,1) double {mustBePositive, mustBeInteger} = 1;
    
    lambda (1,1) double {mustBeNonnegative} = 0;

    mode (1,:) char {mustBeMember(mode, {'identity', 'gauss', 'circulargauss', 'picard'} )} = 'identity';

    options.W (:, :) double = [];

    options.gauss_mu (1,1) double = 0;
    options.gauss_std (1,1) double {mustBePositive} = 100;

    options.picard_log_transform (1,1) logical = true;

    options.enable_disp (1,1) logical = false;
    options.disp_picard_idx (1,1) {mustBeNonnegative, mustBeInteger} = 1;
    options.Fs (1,1) double {mustBePositive} = 48000;
    options.dB_lim (1,2) double = [-inf, inf];
    options.deg_lim (1,2) double = [-180, 180];

end

[N, M] = size(X);
assert(numel(w) == N, 'H, w size mismatch N');

% Compute Discrete Fourier transform matrix
ndx = 0:(num_taps - 1);
F = exp( -1i * (w(:) * ndx) ); %[N x num_taps]

% Apply optional weighting matrix
X_in = X;
if ~isempty(options.W)
    if numel(options.W) == N

        assert(all(options.W(:) > 0), 'options.W must be positive');
        W_sqrt = sqrt(options.W(:));
        F = bsxfun(@times, F, W_sqrt);
        X = bsxfun(@times, X, W_sqrt);

    elseif isequal(size(options.W),  [N, N])        

        [W_sqrt, flag] = chol(options.W);
        assert(flag == 0, 'options.W is not positive definite');
        F = W_sqrt * F;
        X = W_sqrt * X;

    else
        error('options.W size unsupported');
    end
end

% Compute regularization matrix Q
if any( strcmp(mode, {'identity', 'gauss', 'circulargauss'} ) )

    if strcmp(mode, 'identity')
    
        q = ones(1, num_taps);
        Q = lambda * diag(q);
    
    elseif strcmp(mode, 'gauss')
    
        q = 1 - exp( - (ndx - options.gauss_mu).^2 / (2 * options.gauss_std.^2) );
        Q = lambda * diag(q);

    elseif strcmp(mode, 'circulargauss')
   
        exp_sum = zeros(1, num_taps);
        z_std3 = ceil(3 * options.gauss_std / num_taps); % 3 standard deviations
        for i = -z_std3:z_std3
            exp_sum = exp_sum + exp( - (ndx - options.gauss_mu + i * num_taps).^2 / (2 * options.gauss_std.^2) );
        end
        q = max(exp_sum) - exp_sum;
        Q = lambda * diag(q);

    else
        error('Unsupported mode');
    end

    % Compute real-valued least squares solution 
    h = real(F' * F + Q) \ real(F' * X); % LHS is Hermitian

elseif strcmp(mode, 'picard')

    [U, S, V] = svd(F, 'econ'); % Descending singular values    
    s_list = diag(S);
    s_list_ascend = flipud(s_list);
    U_ascend = fliplr(U);        
    V_ascend = fliplr(V);
     
    cdf_s = cumsum(s_list_ascend) / sum(s_list_ascend);

    h = zeros([num_taps, M]);
    for m = 1:M
        U_ascend_X_m = abs(U_ascend' * X(:, m)); %[num_taps x 1]
    
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
    
                % Normalize
                % A_max_abs = max(abs(A), [], 2);
                % A = bsxfun(@rdivide, A, A_max_abs);
                % b = b ./ A_max_abs;
    
                % Solve
                options_linprog = optimoptions("linprog", MaxIterations=1000, ConstraintTolerance=1e-10, OptimalityTolerance=1e-10);
                [d_PC, fval, exitflag] = linprog(ones(N_PC, 1), A, b, [], [], zeros(N_PC, 1), inf(N_PC, 1), options_linprog);
            
            end

            d = zeros(size(V, 2), 1);
            d(1:N_PC) = d_PC;
            Q = V_ascend * diag(d) * V_ascend';
            
            if num_taps > N % Augment with smallest eigenvalue paired with null-space
                V_null = null(V');
                Q = Q + V_null * min(s_list_ascend.^2 + d) * V_null'; 
            end

            % Check eigenvalues if num_taps <= N
            % [V_tmp, D_tmp] = eig(F'*F + Q);
            % norm(sort(diag(D_tmp), 'descend') - sort(s_list_ascend.^2 + d, 'descend'))
            ;
    
        else
            Q = zeros(num_taps);
        end
    
        % Compute real-valued least squares solution
        h(:, m) = real(F' * F + Q) \ real(F' * X); % LHS is Hermitian
    end

else
    error('Unsupported mode');
end

% Compute error
err = norm(F * h - X);
%err_obj = norm(F * h - X)^2 + h'*Q*h;

% Plotting
h_figs = [];
if options.enable_disp && coder.target("MATLAB")

    fontsize = 18;
    hz = options.Fs * w / (2 * pi);
    hz_h = logspace(log10(20), log10(options.Fs / 2), 4096 + 1 )';
    H = freqz(h, 1, hz_h, options.Fs);
    H_in = freqz(h, 1, hz, options.Fs);

    h_figs{1} = figure;
    h_figs{1}.Position = [100, 100, [1200, 900] * (3/4) ];
    tiledlayout(3, 1);
    
    % Magnitude
    nexttile;
    X_in_dB = mag2db(abs(X_in));
    H_dB = mag2db(abs(H));
    H_in_dB = mag2db(abs(H_in));
    semilogx(hz, X_in_dB, 'b-', hz_h, H_dB, 'r-', hz, H_in_dB - X_in_dB, 'g-', 'linewidth', 1.5);
    grid on; axis tight; ylim(options.dB_lim);
    set(gca, 'fontsize', fontsize - 1);
    xlabel('Frequency (Hz)', 'fontsize', fontsize);
    ylabel('Magnitude (dB)', 'fontsize', fontsize);
    title('Magnitude Response', 'fontsize', fontsize + 1);
    h_lg = legend('Target', 'Filter', 'Residual', 'location', 'northwest', 'NumColumns', 3);
    set(h_lg, 'fontsize', fontsize - 1);

    % Phase
    nexttile;
    X_in_deg = rad2deg(angle(X_in));
    H_deg = rad2deg(angle(H));
    H_in_deg = rad2deg(angle(H_in));
    semilogx(hz, X_in_deg, 'b-', hz_h, H_deg, 'r-', hz, H_in_deg - X_in_deg, 'g-', 'linewidth', 1.5);
    grid on; axis tight; ylim(options.deg_lim);
    set(gca, 'fontsize', fontsize - 1);
    xlabel('Frequency (Hz)', 'fontsize', fontsize);
    ylabel('Phase (Degree)', 'fontsize', fontsize);
    title('Phase Response', 'fontsize', fontsize + 1);
    h_lg = legend('Target', 'Filter', 'Residual', 'location', 'northwest', 'NumColumns', 3);
    set(h_lg, 'fontsize', fontsize - 1);

    % Fitted filter time-domain
    nexttile;
    yyaxis left;
    plot(ndx, h, '-', 'linewidth', 1.5);
    ylabel('Filter Taps', 'fontsize', fontsize);
    grid on; axis tight;

    yyaxis right;
    plot(ndx, real(diag(Q)), '-', 'linewidth', 1.5);    
    xlabel('Time Samples', 'fontsize', fontsize);
    ylabel('Regularization Variance', 'fontsize', fontsize);
    title('Time Domain Response', 'fontsize', fontsize + 1);
    set(gca, 'fontsize', fontsize - 1);
    grid on; axis tight;

    if strcmp(mode, 'picard')

        % Picard plot 
        [U, S, V] = svd(F, 'econ'); % Descending singular values
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
            'linewidth', 1.5, 'MarkerSize', 4);
        grid on; axis tight;
        xlabel('Descending Singular Value Index i', 'fontsize', fontsize);
        ylabel('Magnitude', 'fontsize', fontsize);
        title('Picard Plot', 'fontsize', fontsize + 1);
        set(gca, 'fontsize', fontsize - 1);
        if M == 1
            h_lg = legend('$\sigma_i$', '$|u_i^H x|$', '$|u_i^H x| / \sigma_i$',  ...
                '$(\sigma_i^2 + d_i) / \sigma_i$', ...
                '$\sqrt{d_i}$', ...
                '$|u_i^H x| \frac{ \sigma_i} {\sigma_i^2 + d_i}$', ...
                'location' ,'best', 'interpreter', 'latex', 'NumColumns', 2 ); 
        end
        set(h_lg, 'fontsize', fontsize + 1);
    
        nexttile;
        semilogy(1:N_s_list, cdf_s, 'r*-', 1:N_s_list, cdf_ux, 'bs-',  ...
            1:N_s_list, cdf_sr, 'md-', 'linewidth', 1.5, 'MarkerSize', 4); hold on;
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
