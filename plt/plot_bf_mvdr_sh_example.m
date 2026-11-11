function plot_bf_mvdr_sh_example(mode, N, kappa, options)
% Plot examples of bf_mvdr_sh.m

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%mode:          String, array type {'uni_circ', 'spis_circ'}
%N:             Number of array elements
%kappa:         Max condition number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare omni-directivity transducer to spherical piston

%plot_bf_mvdr_sh_example('uni_circ', 2);
%plot_bf_mvdr_sh_example('spis_circ', 2);

%plot_bf_mvdr_sh_example('uni_line', 3);
%plot_bf_mvdr_sh_example('uni_circ', 3);
%plot_bf_mvdr_sh_example('spis_circ', 3);

%plot_bf_mvdr_sh_example('uni_circ', 4);
%plot_bf_mvdr_sh_example('spis_circ', 4);

%plot_bf_mvdr_sh_example('uni_circ', 4, 100);
%plot_bf_mvdr_sh_example('spis_circ', 4, 100);

arguments
    mode (1,:) char {mustBeMember(mode, {'uni_line', 'uni_circ', 'spis_circ'})} = 'uni_circ';
    N (1,1) double {mustBePositive, mustBeInteger} = 4;
    kappa (1,1) double {mustBePositive} = inf;

    options.enable_export (1,1) logical = false;
    options.export_name (1,:) char = [mode, '_N_', num2str(N), '_kappa_', num2str(kappa)];
end

Fs = 24000;
M = 512;
radius = 0.1;
delay_dist = 1;
freq = logspace(log10(20), log10(Fs), M + 1); freq = freq(1:end-1);

spis_arc_len = 0.0254;

theta = pi/2;
phi = 0;

mode_mvdr = 'uniform';

mode_kappa = 'diag';
%mode_kappa = 'thres';
%mode_kappa = 'exp';
%mode_kappa = 'sse';
%mode_kappa = 'ase';

if any(strcmp(mode, {'uni_line', 'uni_circ'}))

    [C_arr, center, h_arr_figs] = sh_enc_arr(0, mode, N, radius, delay_dist, freq, ...
        'noncentered', true, 'enable_disp', true, 'disp_dB_lim', [-3, 3]);
    wavnum = 2 * pi * freq / 343;

elseif strcmp(mode, 'spis_circ')

    max_odr = 50;
    [C_arr, center, h_arr_figs] = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, ...
        'spis_arc_len', spis_arc_len, 'enable_disp', true, 'disp_dB_lim', [-36, 12]);
    wavnum = [];

else
    error('Unsupported mode');
end

[W, fval, SINR, G, h_bf_figs] = bf_mvdr_sh(C_arr, false, freq, theta, phi, kappa, mode_mvdr, ...
    'mode_kappa', mode_kappa, 'center', center, 'wavnum', wavnum, 'enable_disp', true, 'disp_dB_lim', [-48, 14]);

h_figs = [h_arr_figs, h_bf_figs];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.enable_export
    out_dir = 'figs/figs_bf/mvdr';
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    
    for n = 1:numel(h_figs)
        exportgraphics(h_figs{n}, fullfile(out_dir,  [options.export_name, '_', num2str(n), '.png'] ));
    end
end