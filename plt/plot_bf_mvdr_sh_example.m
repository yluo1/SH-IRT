function plot_bf_mvdr_sh_example(mode, N, kappa, mode_kappa, options)
% Plot examples of bf_mvdr_sh.m

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%mode:          String, array type {'uni_line', 'uni_circ', 'spis_circ', 'spis_custom'}
%N:             Number of array elements
%kappa:         Max condition number
%mode_kappa:    String, condition number constraint method {'diag', 'thres', 'exp', 'sse', 'ase', 'minimax'}

%options:       struct

%options.max_arr_gain_pow_dB:      Max array gain bounds in dB power
%options.enable_export:            Logical, if true, export figure

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

% enable_export = true;
% plot_bf_mvdr_sh_example('uni_circ', 5, inf, 'diag', 'max_arr_gain_pow_dB', 6, 'enable_export', enable_export);
% plot_bf_mvdr_sh_example('uni_circ', 5, inf, 'thres', 'max_arr_gain_pow_dB', 6, 'enable_export', enable_export);
% plot_bf_mvdr_sh_example('uni_circ', 5, inf, 'exp', 'max_arr_gain_pow_dB', 6, 'enable_export', enable_export);
% plot_bf_mvdr_sh_example('uni_circ', 5, inf, 'sse', 'max_arr_gain_pow_dB', 6, 'enable_export', enable_export);
% plot_bf_mvdr_sh_example('uni_circ', 5, inf, 'ase', 'max_arr_gain_pow_dB', 6, 'enable_export', enable_export);
% plot_bf_mvdr_sh_example('uni_circ', 5, inf, 'minimax', 'max_arr_gain_pow_dB', 6, 'enable_export', enable_export);


arguments
    mode (1,:) char {mustBeMember(mode, {'uni_line', 'uni_circ', 'spis_circ', 'spis_custom'})} = 'uni_circ';
    N (1,1) double {mustBePositive, mustBeInteger} = 4;
    kappa (1,1) double {mustBePositive} = inf;
    mode_kappa (1,:) char {mustBeMember(mode_kappa, {'diag', 'thres', 'exp', 'sse', 'ase', 'minimax'} )} = 'diag';

    options.enable_export (1,1) logical = false;
    
    options.max_arr_gain_pow_dB (1,1) double {mustBePositive} = inf;

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

if any(strcmp(mode, {'uni_line', 'uni_circ'}))

    [C_arr, center, h_arr_figs] = sh_enc_arr(0, mode, N, radius, delay_dist, freq, ...
        'noncentered', true, 'enable_disp', true, 'disp_dB_lim', [-3, 3]);
    wavnum = 2 * pi * freq / 343;

elseif strcmp(mode, 'spis_circ')

    max_odr = 50;
    [C_arr, center, h_arr_figs] = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, ...
        'spis_arc_len', spis_arc_len, 'enable_disp', true, 'disp_dB_lim', [-36, 12]);
    wavnum = [];

elseif strcmp(mode, 'spis_custom') % Small array element size

    max_odr = 50;
    [C_arr, center, h_arr_figs] = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, ...
        'spis_arc_len', 0.004, 'custom_theta_phi', [deg2rad(10) * ones(1,N); deg2rad(0:(360/N):(360 - 360/N))], ...
        'enable_disp', true, 'disp_dB_lim', [-36, 12]);
    wavnum = [];

else
    error('Unsupported mode');
end

[W, fval, SINR, G, h_bf_figs] = bf_mvdr_sh(C_arr, false, freq, theta, phi, kappa, mode_mvdr, ...
    'mode_kappa', mode_kappa, 'max_arr_gain', db2pow(options.max_arr_gain_pow_dB), ...
    'center', center, 'wavnum', wavnum, 'enable_disp', true, 'disp_dB_lim', [-48, 14]);

h_figs = [h_arr_figs, h_bf_figs];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.enable_export
    out_dir = 'figs/figs_bf/mvdr';
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    export_name = [mode, '_N_', num2str(N), '_', mode_kappa, '_kappa_', num2str(kappa)];

    if isfinite(options.max_arr_gain_pow_dB)   % Augment with finite array gain
        export_name = [export_name, '_arrgain_', num2str(options.max_arr_gain_pow_dB)];
    end
    
    for n = 1:numel(h_figs)
        f_name = [export_name, '_', num2str(n)];
        exportgraphics(h_figs{n}, fullfile(out_dir, [f_name, '.png'] ));
    end
end