function plot_bf_grq_sh_example(mode_sol, tau_tgt_dB, enable_wt, options)
% Plot examples of bf_grq_sh.m

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%mode_sol:      String, solver method   {'spec', 'spec-nowt', 'proj', 'proj-nowt', 'sensitivity', 'sensitivity-nowt', 'efficiency', 'efficiency-nowt'} 
%tau_tgt_dB:    Constant GRQ target (dB)
%enable_wt:     Logical, if true, apply channel weighting

%options:       struct

%options.export_arr:    Logical, if true, export figure 
%options.export_pdf:    Logical, if true, export density
%options.export_tau:    Logical, if true, export tau target
%options.export_bf:     Logical, if true, export beampattern and statistics
%options.export_wt:     Logical, if true, export frequency weights

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare GRQ for varying modes and directivity indices on spherical piston array

%plot_bf_grq_sh_example('sensitivity', db2pow(inf), false);

% For export
% plot_bf_grq_sh_example('sensitivity', inf,  false, 'export_arr', true, 'export_pdf', true, 'export_tau', true, 'export_bf', true);
% plot_bf_grq_sh_example('sensitivity', 6,    false, 'export_tau', true);
% plot_bf_grq_sh_example('sensitivity', 0,    false, 'export_tau', true);
% plot_bf_grq_sh_example('sensitivity', -12,  false, 'export_tau', true);

% plot_bf_grq_sh_example('spec',        6, false, 'export_bf', true);
% plot_bf_grq_sh_example('proj',        6, false, 'export_bf', true);
% plot_bf_grq_sh_example('sensitivity', 6, false, 'export_bf', true);
% plot_bf_grq_sh_example('efficiency',  6, false, 'export_bf', true);

% plot_bf_grq_sh_example('spec',        0, false, 'export_bf', true);
% plot_bf_grq_sh_example('proj',        0, false, 'export_bf', true);
% plot_bf_grq_sh_example('sensitivity', 0, false, 'export_bf', true);
% plot_bf_grq_sh_example('efficiency',  0, false, 'export_bf', true);

% Frequency weighted
% plot_bf_grq_sh_example('sensitivity', inf,  true, 'export_arr', true, 'export_pdf', true, 'export_tau', true, 'export_bf', true, 'export_wt', true);
% plot_bf_grq_sh_example('sensitivity', 6,    true, 'export_tau', true);
% plot_bf_grq_sh_example('sensitivity', 0,    true, 'export_tau', true);
% plot_bf_grq_sh_example('sensitivity', -12,  true, 'export_tau', true);

% plot_bf_grq_sh_example('spec',        6, true, 'export_bf', true);
% plot_bf_grq_sh_example('proj',        6, true, 'export_bf', true);
% plot_bf_grq_sh_example('sensitivity', 6, true, 'export_bf', true);
% plot_bf_grq_sh_example('efficiency',  6, true, 'export_bf', true);

% plot_bf_grq_sh_example('spec',        0, true, 'export_bf', true);
% plot_bf_grq_sh_example('proj',        0, true, 'export_bf', true);
% plot_bf_grq_sh_example('sensitivity', 0, true, 'export_bf', true);
% plot_bf_grq_sh_example('efficiency',  0, true, 'export_bf', true);

arguments 
    mode_sol (1,:) char {mustBeMember(mode_sol, {'spec', 'spec-nowt', 'proj', 'proj-nowt', 'sensitivity', 'sensitivity-nowt', 'efficiency', 'efficiency-nowt'} )} = 'sensitivity';
    tau_tgt_dB (1,:) double = db2pow(6);
    enable_wt (1,1) logical = false;

    options.export_arr (1,1) logical = false;
    options.export_pdf (1,1) logical = false;
    options.export_tau (1,1) logical = false;
    options.export_bf (1,1) logical = false;
    options.export_wt (1,1) logical = false;

end

Fs = 24000;
M = 512;
N = 4;
radius = 0.0254 * 3; % 3 inch
delay_dist = 1;
freq = logspace(log10(20), log10(Fs), M + 1); freq = freq(1:end-1);
is_real = false;

spis_arc_len = 0.0254;
max_odr = 30;

% Generate array
[C_arr, ~, h_arr_figs] = sh_enc_arr(max_odr, 'spis_circ', N, radius, delay_dist, freq, ...
    'spis_arc_len', spis_arc_len, 'enable_disp', true, 'disp_dB_lim', [-36, 12]);

% Couple transducers 2, 4
C_arr = cat(2, C_arr(:,1,:), C_arr(:,2,:)+C_arr(:,4,:),  C_arr(:,3,:));

% Generate accept density for spatial covariance A, R
[A_pdf, h_A_pdf_fig] = sh_pdf_preset('FwdCenter', 'SqExp', 8, is_real, 'ell', 0.5, 'enable_disp', true);
[R_pdf, h_R_pdf_fig] = sh_pdf_preset('FwdCenter', 'SqExpLeftRight', 8, is_real, 'ell', 0.75, 'enable_disp', true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute beamformers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if enable_wt

    % Frequency weights
    [wt, h_fig_wt] = ft_freq_wt_chs(freq, [4000 8000; 100, 10000; 0 1000], [24 0; 12 6; 24 24], 'enable_disp', true, 'disp_dB_lim', [-48, 6], 'soft_knee', true, 'norm_sum', false);
    varargin = {'A_pdf', A_pdf, 'R_pdf', R_pdf, 'enable_disp', true, 'disp_dB_lim', [-48, inf], 'disp_tau_dB_lim', [-36, inf], 'disp_DI_dB_lim', [-inf, inf], 'disp_EP_dB_lim', [-48, inf], 'disp_spin', true};
    [~, f_fig_bf] = bf_grq_sh(C_arr, false, freq, db2pow(tau_tgt_dB), mode_sol, varargin{:}, 'wt', wt);

else

    varargin = {'A_pdf', A_pdf, 'R_pdf', R_pdf, 'enable_disp', true, 'disp_dB_lim', [-48, inf], 'disp_tau_dB_lim', [-inf, inf], 'disp_DI_dB_lim', [-inf, inf], 'disp_EP_dB_lim', [-48, inf], 'disp_spin', true};
    [~, f_fig_bf] = bf_grq_sh(C_arr, false, freq, db2pow(tau_tgt_dB), mode_sol, varargin{:});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.export_arr || options.export_pdf || options.export_bf || options.export_tau || options.export_wt
    
    out_dir = 'figs/figs_bf/grq';
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    f_name = ['bf_', mode_sol, '_tau_', num2str(tau_tgt_dB)];
    if enable_wt
        f_name = [f_name, '_wt'];
    end

    if options.export_arr
        for n = 1:numel(h_arr_figs)
            f_name_arr = ['arr_', num2str(n)];
            exportgraphics(h_arr_figs{n}, fullfile(out_dir, [f_name_arr, '.png'] ));
        end
    end

    if options.export_pdf
        exportgraphics(h_A_pdf_fig, fullfile(out_dir, ['A_pdf', '.png'] ));
        exportgraphics(h_R_pdf_fig, fullfile(out_dir, ['R_pdf', '.png'] ));
    end
    
    if options.export_tau
         exportgraphics(f_fig_bf{1}, fullfile(out_dir, [f_name, '_tau_tgt', '.png'] ));            
    end

    if options.export_bf
        exportgraphics(f_fig_bf{2}, fullfile(out_dir, [f_name, '_bp', '.png'] ));            
        exportgraphics(f_fig_bf{3}, fullfile(out_dir, [f_name, '_bs', '.png'] ));            
        exportgraphics(f_fig_bf{4}, fullfile(out_dir, [f_name, '_spin', '.png'] ));
    end

    if options.export_wt && enable_wt
        exportgraphics(h_fig_wt, fullfile(out_dir, ['freq_wt', '.png'] ));
    end


end