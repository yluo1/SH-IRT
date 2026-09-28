function plot_sh_fit_svd_example
%Plot sh_fit_svd modes

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate and sample from random field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(441);
max_odr = 6;
is_real = false;
N_C = (max_odr + 1)^2;
N_pts = (max_odr + 1)^2;

%Reference field
C_ref = sh_rand(max_odr, 1, is_real);

%Sample at spherical coordinates
%[theta, phi] = sh_grd_fib(N_pts);       %Fibonnaci
[theta, phi] = sh_grd_rand(N_pts);  %Random over sphere

X = sh_dec(C_ref, theta, phi, is_real);
%Add noise
X = X + (randn(size(X)) + randn(size(X)) * 1i) * 5e-2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit across various modes and truncation fractions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mode_list = {'max', 'totalvar', 'bottom'};
trunc_frac_list = [0.01, 0.1, 0.2];

N_trunc = numel(trunc_frac_list);
N_modes = numel(mode_list);

C_list = zeros([N_C, N_modes, N_trunc]);
err_list = zeros(N_modes, N_trunc);
trunc_percent_list = zeros(N_modes, N_trunc);

for i = 1:N_modes
    for j = 1:N_trunc
        [C_list(:, i, j), ~, trunc_percent_list(i, j)] = sh_fit_svd(X, theta, phi, max_odr, is_real, trunc_frac_list(j), mode_list{i});
        err_list(i, j) = err_NSHMSQ(C_ref, C_list(:, i, j)); %Normalized error
    end
end

% Non-regularized least squares
[C_ls]  = sh_fit_svd(X, theta, phi, max_odr, is_real, 0, 'max', true);

%Picard cross criterion
trunc_frac_pc = 0.5;
[C_pc, ~, ~, h_pc]  = sh_fit_svd(X, theta, phi, max_odr, is_real, trunc_frac_pc, 'picardcross', true);
h_picard = h_pc{2};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit after truncating smallest singular values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trunc_frac_all = (0:N_C) ./ N_C;
N_frac_all = N_C + 1;
err_frac_all = zeros(N_frac_all, 1);
trunc_percent_frac_all = zeros(N_frac_all, 1);
for i = 1:N_frac_all
    [C_i, ~, trunc_percent_frac_all(i)] = sh_fit_svd(X, theta, phi, max_odr, is_real, trunc_frac_all(i), 'bottom');
    err_frac_all(i) = err_NSHMSQ(C_ref, C_i); %Normalized error
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot normalized error
fontsize = 16;
h_err = figure;
h_err.Position = [100, 100, 800, 400];
semilogy(100 * trunc_frac_all, err_frac_all, '*-', 'linewidth', 1.5); grid on; axis tight;
xlabel('Bottom Percent Singular Values Removed', 'fontsize', fontsize);
ylabel('Normalized Error', 'fontsize', fontsize);
title('Normalized Fitting Error for Varying Number of Singular Values', 'fontsize', fontsize + 1);
set(gca, 'fontsize', fontsize - 1);

%Plot reference
dB_lim = [-32, 24];
h_ref = sh_plt(C_ref, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Reference', 'disp_theta_phi', [theta, phi]);

%least squares
h_ls = sh_plt(C_ls, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Non-Regularized Least Squares', 'disp_theta_phi', [theta, phi]);

%Picard cross
h_pc = sh_plt(C_pc, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', ['Picard Cross \tau = ', num2str(trunc_frac_pc)], 'disp_theta_phi', [theta, phi]);

%Plot fits
h_list = cell([N_modes, N_trunc]);
for i = 1:N_modes
    for j = 1:N_trunc
        h_ij = sh_plt(C_list(:, i, j), 'mercator', is_real, 'dB_lim', dB_lim, ...
                'title_name', ['TSVD ', mode_list{i}, ' \tau = ', num2str(trunc_frac_list(j)), ', ', num2str(trunc_percent_list(i, j)), '% Cut' ]);
        h_list{i, j} = h_ij{1};
    end
end

err_list

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_dir = 'figs/figs_sh/fit';
if ~isfolder(out_dir)
    mkdir(out_dir);
end

% exportgraphics(h_err, fullfile(out_dir, 'svd_err.png'));
% exportgraphics(h_ref{1}, fullfile(out_dir, 'svd_ref.png'));
% exportgraphics(h_ls{1}, fullfile(out_dir, 'svd_ls.png'));
% exportgraphics(h_pc{1}, fullfile(out_dir, 'svd_pc.png'));
% exportgraphics(h_picard, fullfile(out_dir, 'svd_picard.png'));
% 
% for i = 1:N_modes
%     for j = 1:N_trunc
%         exportgraphics(h_list{i, j}, fullfile(out_dir, ['svd_fit_', num2str(i), '_', num2str(j), '.png']));
%     end
% end
% 
