function tst_sh_pdf_transport(mode, preset_0_dir_name, preset_0_func_name, preset_1_dir_name, preset_1_func_name, t, options)
%Test sh_pdf_transport.m

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%mode:      String, test case, {'EarthMoverDist', 'SlicedWass'}
%preset_0:  String, name of preset, see sh_pdf_preset.m
%preset_1:  String, name of preset, see sh_pdf_preset.m
%t:         [1 x T] Normalized time between 0 and 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Transform between presets

%t =  [0.25, 0.5, 0.75];
%enable_export = true;

% %Presets
%tst_sh_pdf_transport('SlicedWass', 'FwdLeft', 'SqExp', 'FwdRight', 'SqExp',  t, 'preset_title', 'Left to Right', 'enable_export', enable_export);
%tst_sh_pdf_transport('SlicedWass', 'FwdTop',  'SqExp', 'FwdBot',  'SqExp', t, 'preset_title', 'Top to Bottom', 'enable_export', enable_export);
%tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'SqExp', 'FwdCenter', 'SqExpLeftRight',  t, 'preset_title', 'Center to Left + Right', 'enable_export', enable_export);
%tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'SqExpLeftRight', 'FwdCenter', 'SqExpNorthSouth',  t, 'preset_title', 'Left + Right to Poles', 'enable_export', enable_export);
%tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'SqExp', 'FwdCenter', 'DiracRandom', t, 'preset_title', 'Center to Random', 'enable_export', enable_export);
%tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'DiracRandom', 'FwdCenter', 'DiracRandom',  t, 'rseed', 3561, 'preset_title', 'Random to Random', 'enable_export', enable_export);

% %Other
%tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'SqExpDiagLeft', 'FwdCenter', 'SqExpDiagRight',  t, 'preset_title', 'Diagonal to Diagonal', 'enable_export', enable_export);

%tst_sh_pdf_transport('SlicedWass', 'North', 'SqExp', 'North', 'DiracRandom', t, 'preset_title', 'Top to Random', 'enable_export', enable_export);
%tst_sh_pdf_transport('SlicedWass', 'South', 'SqExp', 'South', 'DiracRandom', t, 'preset_title', 'Bottom to Random', 'enable_export', enable_export);
%tst_sh_pdf_transport('SlicedWass', 'SouthWest', 'SqExp', 'SouthWest', 'DiracRandom', t, 'preset_title', 'Bottom Left to Random', 'enable_export', enable_export);
%tst_sh_pdf_transport('SlicedWass', 'NorthWest', 'SqExp', 'NorthWest', 'DiracRandom', t, 'preset_title', 'Top Left to Random', 'enable_export', enable_export);

%tst_sh_pdf_transport('SlicedWass', 'Behind', 'SqExp', 'FwdCenter', 'DiracRandom', t, 'preset_title', 'Behind to Random', 'enable_export', enable_export);
%tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'SqExpLeftRight', 'FwdCenter', 'DiracRandom', t, 'preset_title', 'Left + Right to Random', 'enable_export', enable_export);
%tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'Uniform', 'FwdCenter', 'SqExp',  t, 'preset_title', 'Uniform to Center', 'enable_export', enable_export);

% %SW fit modes:
%tst_sh_pdf_transport('SlicedWass', 'FwdLeft', 'Dirac', 'FwdRight', 'Dirac',  [0.25, 0.5], ...
% 'max_odr', 4, 'SW_fit_mode', 'LeastSquares', ...
%'disp_phase', true, 'fig_size', [450, 600] * (0.75), ...
%'preset_title', 'Left to Right', 'enable_export', false);

%tst_sh_pdf_transport('SlicedWass', 'FwdLeft', 'Dirac', 'FwdRight', 'Dirac',  [0.25, 0.5], ...
% 'max_odr', 4, 'SW_fit_mode', 'NNLS', ...
%'disp_phase', true, 'fig_size', [450, 600] * (0.75), ...
%'preset_title', 'Left to Right', 'enable_export', false);

% tst_sh_pdf_transport('SlicedWass', 'FwdLeft', 'Dirac', 'FwdRight', 'Dirac',  [0.25, 0.5], ...
% 'max_odr', 4, 'SW_fit_mode', 'LS_MSq', 'SW_msq_mode', 'SOMS', ...
% 'disp_phase', true, 'fig_size', [450, 600] * (0.75), ...
% 'preset_title', 'Left to Right', 'enable_export', false);

% tst_sh_pdf_transport('SlicedWass', 'FwdLeft', 'Dirac', 'FwdRight', 'Dirac',  [0.25, 0.5], ...
% 'max_odr', 4, 'SW_fit_mode', 'NNLS_MSq', 'SW_msq_mode', 'SOMS', ...
% 'disp_phase', true, 'fig_size', [450, 600] * (0.75), ...
% 'preset_title', 'Left to Right', 'enable_export', false);

% %Compare across modes
% t =  [0.25, 0.5, 0.75];
% tst_sh_pdf_transport('LinearInterp', 'FwdLeft', 'SqExp', 'FwdRight', 'SqExp',  t, 'preset_title', 'Left to Right LinearInterp');
% tst_sh_pdf_transport('GeometricInterp', 'FwdLeft', 'SqExp', 'FwdRight', 'SqExp',  t, 'preset_title', 'Left to Right GeometricInterp');
% tst_sh_pdf_transport('EarthMoverDist', 'FwdLeft', 'SqExp', 'FwdRight', 'SqExp',  t, 'preset_title', 'Left to Right EarthMoverDist');
% tst_sh_pdf_transport('SlicedWass', 'FwdLeft', 'SqExp', 'FwdRight', 'SqExp',  t, 'preset_title', 'Left to Right SlicedWass');

% %Compare upsampled densities
% t =  [0.25, 0.5, 0.75];
% tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'SqExp', 'FwdCenter', 'SqExpLeftRight',  t, 'rseed', 3561, 'preset_title', 'Center to Left + Right',    'dB_lim', [-48, -6], 'max_odr', 6);
% tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'SqExp', 'FwdCenter', 'SqExpLeftRight',  t, 'rseed', 3561, 'preset_title', 'Center to Left + Right Up', 'dB_lim', [-48, -6], 'max_odr', 6, 'max_odr_upsample', 8);
% tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'SqExp', 'FwdCenter', 'DiracRandom',  t, 'rseed', 3561, 'preset_title', 'Center to Random',    'dB_lim', [-48, -6], 'max_odr', 6);
% tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'SqExp', 'FwdCenter', 'DiracRandom',  t, 'rseed', 3561, 'preset_title', 'Center to Random Up', 'dB_lim', [-48, -6], 'max_odr', 6, 'max_odr_upsample', 8);

% tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'Dirac', 'FwdCenter', 'SqExpLeftRight',  t, 'rseed', 3561, 'preset_title', 'Dirac to Left + Right',    'dB_lim', [-48, -6], 'max_odr', 6);
% tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'Dirac', 'FwdCenter', 'SqExpLeftRight',  t, 'rseed', 3561, 'preset_title', 'Dirac to Left + Right Up', 'dB_lim', [-48, -6], 'max_odr', 6, 'max_odr_upsample', 8);
% tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'Dirac', 'FwdCenter', 'DiracRandom',  t, 'rseed', 3561, 'preset_title', 'Dirac to Random',    'dB_lim', [-48, -6], 'max_odr', 6);
% tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'Dirac', 'FwdCenter', 'DiracRandom',  t, 'rseed', 3561, 'preset_title', 'Dirac to Random Up', 'dB_lim', [-48, -6], 'max_odr', 6, 'max_odr_upsample', 8);

% Compare oversampling
% t =  [0.25, 0.5, 0.75];
% tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'SqExp', 'FwdCenter', 'DiracRandom',  t, 'rseed', 3561, 'preset_title', 'Center to Random Up NC=12', 'dB_lim', [-48, -6]);
% tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'SqExp', 'FwdCenter', 'DiracRandom',  t, 'rseed', 3561, 'preset_title', 'Center to Random Up NC=1', 'dB_lim', [-48, -6], 'SW_N_C_fac', 1);
% tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'SqExp', 'FwdCenter', 'DiracRandom',  t, 'rseed', 3561, 'preset_title', 'Center to Random Up NU=12', 'dB_lim', [-48, -6], 'SW_N_u', 24);

% Compare sliced projection weightings
% t =  [0.25, 0.5, 0.75];
% tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'SqExp', 'FwdCenter', 'DiracRandom',  t, 'rseed', 3561, 'preset_title', 'Center to Random Uniform Wt', 'dB_lim', [-48, -6]);
% tst_sh_pdf_transport('SlicedWass', 'FwdCenter', 'SqExp', 'FwdCenter', 'DiracRandom',  t, 'rseed', 3561, 'preset_title', 'Center to Random Exp. Wt', 'dB_lim', [-48, -6], 'SW_wt_mode', 'exp');

arguments
    mode (1,:) char {mustBeMember(mode,  {'LinearInterp', 'GeometricInterp', 'EarthMoverDist', 'SlicedWass'})} = 'SlicedWass';
    
    preset_0_dir_name (1,:) char  = 'FwdLeft';
    preset_0_func_name (1,:) char = 'SqExp';

    preset_1_dir_name (1,:) char  = 'FwdRight';
    preset_1_func_name (1,:) char  = 'SqExp';

    t (1,:) double {mustBeNonnegative} = [0.25, 0.5, 0.75];

    options.ell (1,1) double {mustBePositive} = 3/4;

    options.max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 8;
    options.max_odr_upsample (1,1) double {mustBeNonnegative, mustBeInteger} = 8; 
    options.rseed (1,1) double {mustBeNonnegative, mustBeInteger} = 28347;

    options.dB_lim (1,2) double = [-48, 0];
    options.fig_pos (1,2) double = [0, 400];
    options.fig_size (1,2) double {mustBePositive} = [450, 300] * (0.75);
    options.disp_phase (1,1) logical = false;

    options.disp_title (1,1) logical = true;
    options.disp_xlabel_axis (1,1) logical = true;

    options.preset_title (1,:) char = [];

    options.SW_fit_mode (1,:) char {mustBeMember(options.SW_fit_mode, {'LeastSquares', 'LS_MSq', 'NNLS', 'NNLS_MSq'})} = 'NNLS_MSq';
    options.SW_msq_mode (1,:) char {mustBeMember(options.SW_msq_mode, {'MS', 'SOMS'})} = 'SOMS';
    options.SW_N_C_fac (1,1) double {mustBePositive} = 12;
    options.SW_N_u (1,1) double {mustBePositive} = 12;
    options.SW_tol (1,1) double {mustBeNonnegative} = 1e-8;
    options.SW_wt_mode (1,:) char {mustBeMember(options.SW_wt_mode, {'uniform', 'exp'})} = 'uniform';

    options.enable_export (1,1) logical = false;

end

is_real = false;
max_odr = options.max_odr;
max_odr_upsample = options.max_odr_upsample;

%Generate expansions at t = 0, 1
rng(options.rseed);
C_pdf = sh_pdf_preset(preset_0_dir_name, preset_0_func_name, max_odr, is_real, 'ell', options.ell);
D_pdf = sh_pdf_preset(preset_1_dir_name, preset_1_func_name, max_odr, is_real, 'ell', options.ell);

if max_odr ~= max_odr_upsample
    C_pdf = sh_resize(C_pdf, max_odr_upsample);
    D_pdf = sh_resize(D_pdf, max_odr_upsample);
end

%Transport
if strcmp(mode, 'EarthMoverDist')
        
    N_val = 200; %Number of discretization points
    K_pdf = [];
    %K_pdf = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, 0, 0, 0.2, is_real), is_real), 'Sum');
    
    
    [E_pdf] = sh_pdf_transport(C_pdf, D_pdf, t, 'EarthMoverDist', is_real, 'EMD_N_val', N_val, 'EMD_K_pdf', K_pdf);

elseif strcmp(mode, 'SlicedWass')

    [E_pdf] = sh_pdf_transport(C_pdf, D_pdf, t, 'SlicedWass', is_real, ...
        'SW_N_C_fac', options.SW_N_C_fac, 'SW_N_u', options.SW_N_u, 'SW_tol', options.SW_tol, ...
        'SW_fit_mode', options.SW_fit_mode, 'SW_msq_mode', options.SW_msq_mode, ...
        'SW_wt_mode', options.SW_wt_mode);

elseif strcmp(mode, 'LinearInterp')

    [E_pdf] = sh_pdf_transport(C_pdf, D_pdf, t, 'LinearInterp', is_real);


elseif strcmp(mode, 'GeometricInterp')

    [E_pdf] = sh_pdf_transport(C_pdf, D_pdf, t, 'GeometricInterp', is_real);

else
    error('Unsupported mode');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dB_lim = options.dB_lim;
fig_pos = options.fig_pos;
fig_size =  options.fig_size;
disp_phase = options.disp_phase;

disp_title = options.disp_title;
disp_xlabel_axis = options.disp_xlabel_axis;

preset_title = options.preset_title;

h_fig_list = {};

h_fig_list{end+1} = sh_plt(C_pdf, 'mercator', is_real, 'dB_lim', dB_lim, 'disp_phase', disp_phase, 'fig_pos', fig_pos, 'fig_size', fig_size, ...
    'title_name_override', [preset_title, ' (t = 0)'], 'disp_cb', false, 'disp_title', disp_title, ...
    'disp_xlabel', disp_xlabel_axis, 'disp_xticks', disp_xlabel_axis);
fig_pos(1) = fig_pos(1) + fig_size(1);

for i = 1:numel(t)
    h_fig_list{end+1} = sh_plt(E_pdf(:, :, i), 'mercator', is_real, 'dB_lim', dB_lim, 'disp_phase', disp_phase, 'fig_pos', fig_pos,  'fig_size', fig_size, ...
        'title_name_override', [preset_title, ' (t = ', num2str(t(i)), ')'], 'disp_cb', false, 'disp_title', disp_title, ...
        'disp_ylabel', false, 'disp_yticks', false, ...
        'disp_xlabel', disp_xlabel_axis, 'disp_xticks', disp_xlabel_axis);
    
    fig_pos(1) = fig_pos(1) + fig_size(1);
end

h_fig_list{end+1} = sh_plt(D_pdf, 'mercator', is_real, 'dB_lim', dB_lim, 'disp_phase', disp_phase, 'fig_pos', fig_pos,  'fig_size', fig_size, ...
    'title_name_override', [preset_title, ' (t = 1)'], 'disp_cb', false, 'disp_title', disp_title, ...
    'disp_ylabel', false, 'disp_yticks', false, ...
    'disp_xlabel', disp_xlabel_axis, 'disp_xticks', disp_xlabel_axis);
fig_pos(1) = fig_pos(1) + fig_size(1);


%Export figs
if options.enable_export
    out_dir = 'figs_sh';
    mkdir(out_dir);
    for n = 1:numel(h_fig_list)
        export_name = [preset_0_dir_name, preset_0_func_name, '_', preset_1_dir_name, preset_1_func_name, '_' num2str(n)];        
        exportgraphics(h_fig_list{n}{1} , fullfile(out_dir, [export_name, '.png']));
    end
end
;

