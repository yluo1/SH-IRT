function plot_sh_conv_example
%Plot SqExpNorm filters convolved with each other in SH domain

%Author: Yuancheng Luo, 2026

rng(2134);
% P_C = 3;
% P_D = 3;

P_C = 6;
P_D = 6;

theta = pi/2;
phi = 0;
ell = 0.5;
C = sh_enc_rbf('SqExpNorm', P_C, 0, phi, ell, false, false);
%D = sh_rand(P_D, 1);
D = sh_enc_rbf('SqExpNorm', P_D, theta, phi, ell, false, false);
E = sh_conv(D, C);

%Plotting
fig_size = [560, 420 * (3/4)];
varargin = {'dB_lim', [-32, 18], 'disp_phase', false, 'fig_size', fig_size, 'fontsize', 18};

h_C = sh_plt(C, 'mercator', false, 'title_name', 'Spatial Filter', varargin{:}, 'disp_cb', false, 'disp_ylabel', true, 'disp_yticks', true);
h_D = sh_plt(D, 'mercator', false, 'title_name', 'Random', varargin{:}, 'disp_cb', false, 'disp_ylabel', false, 'disp_yticks', false);
h_E = sh_plt(E, 'mercator', false, 'title_name', 'Filtered', varargin{:}, 'disp_cb', true, 'disp_ylabel', false, 'disp_yticks', false);

%Export figures
out_dir = 'figs';
if ~isfolder(out_dir)
    mkdir(out_dir);
end

% exportgraphics(h_C{1}, fullfile(out_dir, 'sample_spatial_SqExP.png'));
% exportgraphics(h_D{1}, fullfile(out_dir, 'sample_spatial_random_func.png'));
% exportgraphics(h_E{1}, fullfile(out_dir, 'sample_spatial_filtered_func.png'));
