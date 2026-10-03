function [C, err, trunc_percent, h_figs] = sh_fit_svd(X, theta, phi, max_odr, is_real, trunc_frac, mode, options)
%Truncated singular value decomposition least squares fit of spherical harmonic bases 
%min_C ||Y(theta, phi) * C - X||^2 for SH bases Y evaluated at theta, phi, coefficients C

%Different modes zero the inverse singular values: inv_singular_values(mask) = 0

%mode = 'max' (recommend trunc_frac = 0.1) truncate singular values below fraction of largest singular value
%mask: singular_values <= max(singular_values) * trunc_frac

%mode = 'totalvar' (recommend trunc_frac =  0.01) truncate singular values below faction of total variance
%mask: cumsum(descending_eigenvalues) / sum(descending_eigenvalues) > 1 - trunc_frac

%mode = 'bottom' (recommend trunc_frac =  0.1) truncate bottom fraction of singular values
%mask: smallest ranked 100 * trunc_frac percent of singular values

%mode = 'picard' (recommend trunc_frac = 0.5) truncate tail singular values when the discrete Picard condition cannot be met
%mask: singular values smaller than 
%first ascending singular value normalized CDF / observation projected on left singular value normalized CDF >= trunc_frac

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X:             [N x M] N measurements of M functions

%theta:         [N x 1] Co-latitude [0, pi]
%phi:           [N x 1] Azimuth [0, 2 * pi)

%max_odr:       Max SH order 
%is_real:       Logical, if true, evaluate real SH

%trunc_frac:    Fraction of largest singular values
%mode:          String, truncation method {'max', 'totalvar', 'bottom', 'picard'}

%options:               struct
%options.W:             [], [N x 1], or [N x N] weighting matrix for weighted least squares
%                       [] uniform weighting
%                       [N x 1] Positive weight vector
%                           min_C ||diag(W.^(1/2)) * (Y(theta, phi) * C - X)||^2
%                       [N x N] Positive definite weight matrix
%                           min_C ||chol(W) * (Y(theta, phi) * C - X)||^2

%options.enable_disp:   Logical, if true, plot fit and Picard plot
%options.disp_dB_lim:   [1 x 2] dB range [min, max]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:                 [(max_odr + 1)^2 x M] SH coefficients
%err:               Scalar, norm( W^(1/2) * (Y * C - X) )
%trunc_percent:     Percentage of singular values truncated
%h_figs:            Handle to figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Least squares SH fit to random function distributed along 
%spherical Fibonnaci points with varying truncation fractions

% rng(441);
% max_odr = 10;
% is_real = false;
% N_pts = (max_odr + 1)^2;
% X = randn(N_pts, 1) + randn(N_pts, 1) * 1i;
% [theta, phi] = sh_grd_fib(N_pts);

% [C_trunc_0, err_trunc_0, trunc_percent_0] = sh_fit_svd(X, theta, phi, max_odr, is_real, 0); err_trunc_0
% sh_plt(C_trunc_0, 'mercator', is_real); trunc_percent_0

% [C_trunc_01, err_trunc_01, trunc_percent_01] = sh_fit_svd(X, theta, phi, max_odr, is_real, 0.01); err_trunc_01
% sh_plt(C_trunc_01, 'mercator', is_real); trunc_percent_01

% [C_trunc_50, err_trunc_50, trunc_percent_50] = sh_fit_svd(X, theta, phi, max_odr, is_real, 0.50); err_trunc_50
% sh_plt(C_trunc_50, 'mercator', is_real); trunc_percent_50

% [C_trunc_90, err_trunc_90, trunc_percent_90] = sh_fit_svd(X, theta, phi, max_odr, is_real, 0.90); err_trunc_90
% sh_plt(C_trunc_90, 'mercator', is_real); trunc_percent_90

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Least squares SH fit to radial basis function along 
%spherical Fibonnaci points with varying max_odr

% rng(441);
% max_odr = 10;
% is_real = false;
% N_pts = (max_odr + 1)^2;
% X = randn(N_pts, 1) + randn(N_pts, 1) * 1i;
% [theta, phi] = sh_grd_fib(N_pts);

% [C_max_odr_9, err_max_odr_9] = sh_fit_svd(X, theta, phi, 9, is_real, 0); 
% sh_plt(C_max_odr_9); err_max_odr_9

% sh_plt(sh_val(max_odr, theta, phi, is_real) \ X, 'mercator', is_real);

% [C_max_odr_11, err_max_odr_11] = sh_fit_svd(X, theta, phi, 11, is_real, 0); 
% sh_plt(C_max_odr_11, 'mercator', is_real); err_max_odr_11

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Weighted least squares SH fit

% rng(441);
% max_odr = 4;
% is_real = false;
% C_ref = sh_rand(max_odr, 1, is_real);
% N_pts = (max_odr + 1)^2;
% [theta, phi] = sh_grd_rand(N_pts);
% X = sh_dec(C_ref, theta, phi, is_real);
% % X = X + (randn(N_pts, 1) + randn(N_pts, 1) * 1i) * 1e-1;

% tau = 0.1;
%C_ls = sh_fit_svd(X, theta, phi, max_odr, is_real, tau, 'max');
%C_diag_W = sh_fit_svd(X, theta, phi, max_odr, is_real, tau, 'max', 'W', rand(N_pts, 1));
%W_psd = rand(N_pts, N_pts); W_psd = W_psd * W_psd';
%C_psd_W = sh_fit_svd(X, theta, phi, max_odr, is_real, tau, 'max', 'W', W_psd);

%sh_plt(C_ref, 'mercator', is_real, 'title_name', 'Reference', 'disp_theta_phi', [theta, phi]);
%sh_plt(C_ls, 'mercator', is_real, 'title_name', 'Uniform Weighting');
%sh_plt(C_diag_W, 'mercator', is_real, 'title_name', 'Diagonal Weighting');
%sh_plt(C_psd_W, 'mercator', is_real, 'title_name', 'PSD Weighting');

arguments
    X (:,:) double {coder.mustBeComplex} = complex(0);

    theta (:,1) double = [0];
    phi   (:,1) double = [0];

    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;

    is_real (1,1) logical = false;

    trunc_frac (1,1) double {mustBeNonnegative, mustBeLessThanOrEqual(trunc_frac, 1)} = 0;

    mode (1,:) char {mustBeMember(mode, {'max', 'totalvar', 'bottom', 'picard'})} = 'max';

    options.W (:, :) double = [];
    options.enable_disp (1,1) logical = false;
    options.disp_dB_lim (1,2) double = [-inf inf];    
end

%Compute SH bases
Y = sh_val(max_odr, theta, phi, is_real);  %[N x (max_odr + 1)^2]
N = size(Y, 1);
M = size(X, 2);

%Apply optional weighting matrix
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

%Compute singular value decomposition
[U, S, V] = svd(Y, "econ");
s_list = diag(S);   %List of singular values
N_s_list = numel(s_list);

%Compute pseudo-inverse singular values
inv_s_list = 1 ./ s_list;

%Compute trunc_mask for setting inv_singular_values(trunc_mask) = 0
if any(strcmp(mode, {'max', 'totalvar', 'bottom'}))

    if strcmp(mode, 'max') 
    
        trunc_mask = s_list <= (max(s_list) * trunc_frac);
    
    elseif strcmp(mode, 'totalvar')
    
        [eig_list_descend, idx] = sort(s_list.^2, 'descend');
        trunc_mask = false([N_s_list, 1]);
        trunc_mask(idx(cumsum(eig_list_descend) / sum(eig_list_descend) > 1 - trunc_frac)) = true;
    
    elseif strcmp(mode, 'bottom')
    
        [s_list_ascend, idx] = sort(s_list, 'ascend');
        trunc_mask = false([N_s_list, 1]);   
        trunc_mask(idx(1:ceil(N_s_list * trunc_frac))) = true;

    else
        error('Unsupported mode');
    end

    inv_s_list(trunc_mask) = 0;
    trunc_percent = sum(trunc_mask) / N_s_list * 100;
    
    inv_S = S;
    inv_S(1:N_s_list, 1:N_s_list) = diag(inv_s_list);
    C = V * inv_S' * U' * X;

elseif strcmp(mode, 'picard')
    
    [s_list_ascend, idx] = sort(s_list, 'ascend');
    U_ascend = U(:, idx);     
    cdf_s = cumsum(s_list_ascend) / sum(s_list_ascend);

    C = zeros([(max_odr + 1)^2, M]);
    trunc_percent_list = zeros(1, M);
    for m = 1:M
        U_ascend_X_m = abs(U_ascend' * X(:, m)); %[N_s_list x 1]
    
        cdf_ux = cumsum(U_ascend_X_m) / sum(U_ascend_X_m);
    
        idx_picard_cross = find((cdf_s ./ cdf_ux) >= trunc_frac, 1, 'first');
        trunc_mask = false([N_s_list, 1]);  
        trunc_mask(1:(idx_picard_cross-1)) = true;
        trunc_mask = flipud(trunc_mask);
    
        inv_s_list(trunc_mask) = 0;
        trunc_percent_list(m) = sum(trunc_mask) / N_s_list * 100;
        
        inv_S = S;
        inv_S(1:N_s_list, 1:N_s_list) = diag(inv_s_list);
        C(:, m) = V * inv_S' * U' * X(:, m);
    end
    trunc_percent = mean(trunc_percent_list);

else
    error('Unsupported mode');
end

%Compute error
err = norm(Y * C - X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_figs = [];
if options.enable_disp && coder.target("MATLAB")

    % Plot fit
    h_figs = sh_plt(C, 'mercator', is_real, 'disp_theta_phi', [theta, phi], 'title_name', ['SVD Fit ', mode, ' \tau = ', num2str(trunc_frac)], 'dB_lim', options.disp_dB_lim);

    % Picard plot
    [~, idx_descend] = sort(s_list, 'descend'); %Descending singular values
    U_descend = U(:, idx_descend);
    s_list_descend = s_list(idx_descend);
    U_descend_X_mu = mean(abs(U_descend' * X), 2); %[N_s_list x 1]

    cdf_s = cumsum(flipud(s_list_descend(:))) / sum(s_list_descend);
    cdf_ux = cumsum(flipud(U_descend_X_mu(:))) / sum(U_descend_X_mu);
    
    idx_picard_cross = find( cdf_s>= cdf_ux, 1, 'first');
    idx_picard_cross_trunc = find((cdf_s ./ cdf_ux) >= trunc_frac, 1, 'first');

    fontsize = 16;    
    h_figs{end+1} = figure;
    tiledlayout(2, 1);
    nexttile;
    semilogy(1:N_s_list, s_list_descend, 'r*-', 1:N_s_list, U_descend_X_mu, 'bs-', ...
        1:N_s_list, U_descend_X_mu ./ s_list_descend, 'o-', 'linewidth', 1.5);
    grid on; axis tight;
    xlabel('Descending Singular Value Index i', 'fontsize', fontsize);
    ylabel('Magnitude', 'fontsize', fontsize);
    title('Picard Plot', 'fontsize', fontsize + 1);
    set(gca, 'fontsize', fontsize - 1);
    if M == 1
        h_lg = legend('$\sigma_i$', '$|u_i^H x|$',       '$|u_i^H x| / \sigma_i$',          'location' ,'best', 'interpreter', 'latex'); 
    else
        h_lg = legend('$\sigma_i$', 'mean($|u_i^H x|$)', 'mean($|u_i^H x|$) $ / \sigma_i$', 'location' ,'best', 'interpreter', 'latex'); 
    end
    set(h_lg, 'fontsize', fontsize - 1);

    nexttile;
    semilogy(1:N_s_list, cdf_s, 'r*-', 1:N_s_list, cdf_ux, 'bs-',  'linewidth', 1.5); hold on;
    xline(idx_picard_cross, 'm--', 'linewidth', 2);
    xline(idx_picard_cross_trunc, 'k--', 'linewidth', 2);
   
    grid on; axis tight;
    xlabel('Ascending Singular Value Index i', 'fontsize', fontsize);
    ylabel('Magnitude', 'fontsize', fontsize);
    title('Normalized Cumulative Distribution Function', 'fontsize', fontsize + 1);
    set(gca, 'fontsize', fontsize - 1);
    if M == 1
        h_lg = legend('CDF($\sigma_i$)', 'CDF($|u_i^H x|$)',       'Picard Cross: $\tau = 1$',  ['Picard Cross: $\tau = ', num2str(trunc_frac), '$'], 'location' ,'best', 'interpreter', 'latex'); 
    else
        h_lg = legend('CDF($\sigma_i$)', 'CDF(mean($|u_i^H x|$))', 'Picard Cross: $\tau = 1$',  ['Picard Cross: $\tau = ', num2str(trunc_frac), '$'], 'location' ,'best', 'interpreter', 'latex'); 
    end
    set(h_lg, 'fontsize', fontsize - 1);

end