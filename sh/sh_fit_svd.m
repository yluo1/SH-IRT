function [C, err, trunc_percent, h_figs] = sh_fit_svd(X, theta, phi, max_odr, is_real, trunc_frac, mode, enable_disp)
%Least-squares spherical harmonic basis fit min ||Y*C - X||^2 
%via truncated singular value decomposition: inv_singular_values(mask) = 0

%mode = 'max' (recommend trunc_frac = 0.1)
%mask: singular_values <= max(singular_values) * trunc_frac

%mode = 'totalvar' (recommend trunc_frac =  0.01)
%mask: cumsum(descending_eigenvalues) / sum(descending_eigenvalues) > 1 - trunc_frac

%mode = 'bottom' (recommend trunc_frac =  0.1)
%mask: smallest ranked 100 * trunc_frac percent of singular values

%mode = 'picardcross' (recommend trunc_frac = 0.5) 
%mask: singular values smaller than 
%first ascending singular value CDF / observation projected on left singular value CDF >= trunc_frac

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X:             [N x M] N measurements of M functions

%theta:         [N x 1] Co-latitude [0, pi]
%phi:           [N x 1] Azimuth [0, 2 * pi)

%max_odr:       Max SH order 
%is_real:       Logical, if true, evaluate real SH

%trunc_frac:    Fraction of largest singular values
%mode:          String, truncation method {'max', 'totalvar', 'bottom', 'picardcross'}
%enable_disp:   Logical, if true, plot fit and Picard plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:                 [(max_odr + 1)^2 x M] SH coefficients
%err:               Scalar, norm(Y*C - X);
%trunc_percent:     Percentage of singular values truncated
%h_figs:            Handle to figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_fit_svd', '-o', 'sh/sh_fit_svd_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Least-squares SH fit to random function distributed along 
%spherical Fibonnaci points with varying truncation fractions

% rng(441);
% max_odr = 10;
% is_real = false;
% N_pts = (max_odr + 1)^2;
% X = randn(N_pts, 1) + randn(N_pts, 1) * 1i;
% [theta, phi] = sh_fib(N_pts);

% [C_trunc_0, err_trunc_0, trunc_percent_0] = sh_fit_svd(X, theta, phi, max_odr, is_real, 0); err_trunc_0
% sh_plt(C_trunc_0); trunc_percent_0

% [C_trunc_01, err_trunc_01, trunc_percent_01] = sh_fit_svd(X, theta, phi, max_odr, is_real, 0.01); err_trunc_01
% sh_plt(C_trunc_01); trunc_percent_01

% [C_trunc_50, err_trunc_50, trunc_percent_50] = sh_fit_svd(X, theta, phi, max_odr, is_real, 0.50); err_trunc_50
% sh_plt(C_trunc_50); trunc_percent_50

% [C_trunc_90, err_trunc_90, trunc_percent_90] = sh_fit_svd(X, theta, phi, max_odr, is_real, 0.90); err_trunc_90
% sh_plt(C_trunc_90); trunc_percent_90

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Least-squares SH fit to radial basis function along 
%spherical Fibonnaci points with varying max_odr

% rng(441);
% max_odr = 10;
% is_real = false;
% N_pts = (max_odr + 1)^2;
% X = randn(N_pts, 1) + randn(N_pts, 1) * 1i;
% [theta, phi] = sh_fib(N_pts);

% [C_max_odr_9, err_max_odr_9] = sh_fit_svd(X, theta, phi, 9, is_real, 0); 
% sh_plt(C_max_odr_9); err_max_odr_9

% sh_plt(sh_val(max_odr, theta, phi, is_real) \ X);

% [C_max_odr_11, err_max_odr_11] = sh_fit_svd(X, theta, phi, 11, is_real, 0); 
% sh_plt(C_max_odr_11); err_max_odr_11


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare Matlab versus Mex

% rng(441);
% max_odr = 10;
% is_real = false;
% N_pts = (max_odr + 1)^2;
% X = randn(N_pts, 1) + randn(N_pts, 1) * 1i;
% [theta, phi] = sh_fib(N_pts);

% C_trunc_mat = sh_fit_svd(X, theta, phi, max_odr, is_real, 0.90); 
% C_trunc_mex = sh_fit_svd_mex(X, theta, phi, max_odr, is_real, 0.90); 
% err = norm(C_trunc_mat - C_trunc_mex)

arguments
    X (:,:) double {coder.mustBeComplex} = complex(0);

    theta (:,1) double = [0];
    phi   (:,1) double = [0];

    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;

    is_real (1,1) logical = false;

    trunc_frac (1,1) double {mustBeNonnegative, mustBeLessThanOrEqual(trunc_frac, 1)} = 0;

    mode (1,:) char {mustBeMember(mode, {'max', 'totalvar', 'bottom', 'picardcross'})} = 'max';

    enable_disp (1,1) logical = false;
end

%Compute SH bases
Y = sh_val(max_odr, theta, phi, is_real);  %[N x (max_odr + 1)^2]

%Compute singular value decomposition
[U, S, V] = svd(Y, "econ");
s_list = diag(S);   %List of singular values
N_s_list = numel(s_list);

%Compute pseudo-inverse singular values
inv_s_list = 1 ./ s_list;

%Compute trunc_mask for setting inv_singular_values(trunc_mask) = 0
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

elseif strcmp(mode, 'picardcross')

    [s_list_ascend, idx] = sort(s_list, 'ascend');
    U_ascend = U(:, idx);
    U_ascend_X_mu = mean(abs(U_ascend' * X), 2); %[N_s_list x 1]
     
    cdf_s = cumsum(s_list_ascend) / sum(s_list_ascend);
    cdf_ux = cumsum(U_ascend_X_mu) / sum(U_ascend_X_mu);

    idx_picard_cross = find((cdf_s ./ cdf_ux) >= trunc_frac, 1, 'first');
    trunc_mask = false([N_s_list, 1]);  
    trunc_mask(1:(idx_picard_cross-1)) = true;
    trunc_mask = flipud(trunc_mask);

else
    error('Unsupported mode');
end

inv_s_list(trunc_mask) = 0;
trunc_percent = sum(trunc_mask) / N_s_list * 100;

inv_S = S;
inv_S(1:N_s_list, 1:N_s_list) = diag(inv_s_list);
C = V * inv_S' * U' * X;

%Compute error
err = norm(Y * C - X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_figs = [];
if enable_disp && coder.target("MATLAB")

    % Plot fit
    h_figs = sh_plt(C, 'mercator', is_real, 'disp_theta_phi', [theta, phi], 'title_name', ['SVD Fit ', mode, ' \tau = ', num2str(trunc_frac)]);

    % Plot Picard plot
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
    semilogy(1:N_s_list, s_list_descend, 'r*-', 1:N_s_list, U_descend_X_mu, 'bs-', 'linewidth', 1.5);
    grid on; axis tight;
    xlabel('Descending Singular Value Index i', 'fontsize', fontsize);
    ylabel('Magnitude', 'fontsize', fontsize);
    title('Picard Plot', 'fontsize', fontsize + 1);
    set(gca, 'fontsize', fontsize - 1);
    h_lg = legend('$\sigma_i$', '$|u_i^H x_i|$', 'location' ,'best', 'interpreter', 'latex'); 
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
    h_lg = legend('CDF($\sigma_i$)', 'CDF($|u_i^H x_i|$)', 'Picard Cross: $\tau = 1$',  ['Picard Cross: $\tau = ', num2str(trunc_frac), '$'], 'location' ,'best', 'interpreter', 'latex'); 
    set(h_lg, 'fontsize', fontsize - 1);

end