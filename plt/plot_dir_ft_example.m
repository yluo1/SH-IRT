function plot_dir_ft_example
%Plot directional filter, random function, and filtered function

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(2134);
% P_C = 3;
% P_D = 3;

%P_C = 12;

% P_C = 6;
% P_D = 6;

P_C = 7;
P_D = 7;

theta = pi/2;
phi = 0;
ell = 0.5;
D = sh_enc_rbf('SqExp', P_D, theta, phi, ell, false, false);
C = sh_rand(P_C, 1);
[E, A] = sh_mul(D, C); 

%Plotting
%fig_size = [560, 420 * (3/4)];
fig_size = [560, 420 * 0.65];

varargin = {'dB_lim', [-32, 18], 'disp_phase', false, 'fig_size', fig_size, 'fontsize', 20};

h_D = sh_plt(D, 'mercator', false, 'title_name', 'Directional Function', varargin{:}, 'disp_cb', false, 'disp_ylabel', true, 'disp_yticks', true);
h_C = sh_plt(C, 'mercator', false, 'title_name', 'Random Function', varargin{:}, 'disp_cb', false, 'disp_ylabel', false, 'disp_yticks', false);
h_E = sh_plt(E, 'mercator', false, 'title_name', 'Product Function', varargin{:}, 'disp_cb', true, 'disp_ylabel', false, 'disp_yticks', false);

%Export figures
out_dir = 'figs';
if ~isfolder(out_dir)
    mkdir(out_dir);
end

exportgraphics(h_D{1}, fullfile(out_dir, 'sample_SqExP_dir.png'));
exportgraphics(h_C{1}, fullfile(out_dir, 'sample_random_func.png'));
exportgraphics(h_E{1}, fullfile(out_dir, 'sample_filtered_func.png'));
