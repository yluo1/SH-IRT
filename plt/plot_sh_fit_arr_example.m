function plot_sh_fit_arr_example
% Plot figures for sample fitted piston SH expansion and pre-processing

% Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dat = load('dat/sample_exp_centered.mat');

is_real = false;
N_C = size(dat.C, 1);
P = sqrt(N_C) - 1;
dB_lim = [-40, 20];

% Plot SH expansion centered on piston surface
h_figs_cp = sh_plt(dat.C(:, 3:end), 'horizontal', is_real, 'hz', dat.freq(3:end), 'dB_lim', dB_lim);

% Translate the expansion center to the baffle centroid
omega = 2 * pi * dat.freq;
[C] = sh_translate(dat.C, omega, is_real, dat.loc, 'max_odr_fit', P + 1);
h_figs_cb = sh_plt(C(:, 3:end), 'horizontal', is_real, 'hz', dat.freq(3:end), 'dB_lim', dB_lim);

% Display cardinal direction responses
% Fs = dat.freq(end) * 2;
% M = size(dat.C, 2);
%t = (0:( 2*(M-1) - 1 ) ) / Fs ;
%sh_plt(dat.C, 'cardinal', false, 'hz', dat.freq);
%sh_plt(dat.C, 'cardinal', false, 'hz', dat.freq, 't', t);

%sh_plt(C, 'cardinal', false, 'hz', dat.freq);
%sh_plt(C, 'cardinal', false, 'hz', dat.freq, 't', t);

% Duplicate and array three expansions rotated in the frontal half-space
delay_dist = 0;
radius = abs(dat.loc(1));
custom_theta_phi = deg2rad([90, -30; 90, 30; 90 0])';
custom_rot_zyx   = deg2rad([-30 0 0; 30 0 0; 0 0 0])';
[C_arr, center, h_figs_arr] = sh_enc_arr(P+1, 'fit_rot_custom', 3, radius, delay_dist, dat.freq, ...
    'custom_C', C,  'custom_theta_phi', custom_theta_phi,  'custom_rot_zyx', custom_rot_zyx, ...
    'enable_disp', true, 'disp_dB_lim', dB_lim);

% Design MVDR beam at 0, 30, 45 degrees
theta_steer = deg2rad(90);
phi_steer_list = deg2rad([0, 30, 45]);
kappa = inf;
mode_mvdr = 'uniform';

bf_mvdr_sh(C_arr(:, :, 3:end), false, dat.freq(3:end), theta_steer, phi_steer_list(1), kappa, mode_mvdr,...
    'enable_disp', true, 'disp_dB_lim', dB_lim);


;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_dir = 'figs/figs_bf/arr';
if ~isfolder(out_dir)
    mkdir(out_dir);
end

% exportgraphics(h_figs_cp{1}, fullfile(out_dir, 'fit_pist_center_pist.png' ));
% exportgraphics(h_figs_cb{1}, fullfile(out_dir, 'fit_pist_center_baffle.png' ));

% exportgraphics(h_figs_arr{1}, fullfile(out_dir, 'fit_pist_arr_N_3_topology.png' ));
% exportgraphics(h_figs_arr{2}, fullfile(out_dir, 'fit_pist_arr_N_3_1.png' ));
% exportgraphics(h_figs_arr{3}, fullfile(out_dir, 'fit_pist_arr_N_3_2.png' ));
% exportgraphics(h_figs_arr{4}, fullfile(out_dir, 'fit_pist_arr_N_3_3.png' ));