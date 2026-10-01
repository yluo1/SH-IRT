function plot_sh_fit_rbf_example
%Plot sh_fit_rbf.m modes

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare RBF fits for varying length-scales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(43);

max_odr = 4;
max_odr_fit = 16;
is_real = false;
N_pts = 50;

% Reference field
C_ref = sh_rand(max_odr, 1, is_real);
C_ref = sh_resize(C_ref, max_odr_fit);

% Sample at spherical coordinates
[theta, phi] = sh_grd_rand(N_pts);  %Random over sphere

X = sh_dec(C_ref, theta, phi, is_real);
% Add noise
noise_std = 5e-1;
X = X + (randn(size(X)) + randn(size(X)) * 1i) * noise_std;

% Plot ref.
dB_lim =  [-60, 20];
h_ref_ell = sh_plt(C_ref, 'mercator', is_real, 'dB_lim', dB_lim, 'disp_theta_phi', [theta, phi], 'disp_theta_phi_markersize', 12);

% Fit RBFs for various length-scales
ell_list = [0.2, 0.4, 0.8];
N_ell = numel(ell_list);
C_fit_list = zeros((max_odr_fit + 1)^2, N_ell);
err_list = zeros(1, N_ell);
h_rbf_ell_list = cell(1, N_ell);
varargin = {'max_iter', 0, 'lb', [0, noise_std^2], 'ub', [1, noise_std^2]};
mode = 'Mat32';
for n = 1:N_ell
    C_fit_list(:, n) = sh_fit_rbf(X, theta, phi, max_odr_fit, is_real, mode, ell_list(n), noise_std^2, varargin{:}); 
    err_list(n) = err_NSHMSQ(C_ref, C_fit_list(:, n));
    h_tmp = sh_plt(C_fit_list(:, n), 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', [mode, ' length-scale = ', num2str(ell_list(n))]);
    h_rbf_ell_list{n} = h_tmp{1};
end
err_list

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare RBF to SVD, Tikhonov regularization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargin = {'max_odr', 5, 'max_odr_fit', 10, 'N_pts', 100};

[err_svd, h_ref, h_svd]     = tst_sh_fit('TSVD', 'NSHMSQ', varargin{:}, 'svd_mode', 'picard', 'svd_trunc_frac', 0.5);
[err_tr, h_ref, h_tr]       = tst_sh_fit('TRegu', 'NSHMSQ', varargin{:}, 'tr_mode', 'picard', 'tr_lambda', 0.5);

[err_SqExp, ~, h_SqExp]     = tst_sh_fit('SqExp', 'NSHMSQ', varargin{:});
[err_Mat52, ~, h_Mat52]     = tst_sh_fit('Mat52', 'NSHMSQ', varargin{:});
[err_Mat32, ~, h_Mat32]     = tst_sh_fit('Mat32', 'NSHMSQ', varargin{:});
[err_Exp, ~, h_Exp]         = tst_sh_fit('Exp', 'NSHMSQ', varargin{:});

err_tr
err_SqExp
err_Mat52
err_Mat32
err_Exp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% out_dir = 'figs/figs_sh/fit';
% if ~isfolder(out_dir)
%     mkdir(out_dir);
% end
% 
% exportgraphics(h_ref_ell{1}, fullfile(out_dir, ['rbf_ell_ref.png']));
% for n = 1:N_ell
%     exportgraphics(h_rbf_ell_list{n}, fullfile(out_dir, ['rbf_ell_', num2str(n), '.png']));
% end

% exportgraphics(h_ref, fullfile(out_dir, 'rbf_ref.png'));
% exportgraphics(h_tr, fullfile(out_dir, 'rbf_tr.png'));
% 
% exportgraphics(h_SqExp, fullfile(out_dir, 'rbf_sqexp.png'));
% exportgraphics(h_Mat52, fullfile(out_dir, 'rbf_mat52.png'));
% exportgraphics(h_Mat32, fullfile(out_dir, 'rbf_mat32.png'));
% exportgraphics(h_Exp, fullfile(out_dir, 'rbf_exp.png'));
