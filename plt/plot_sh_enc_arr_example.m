function plot_sh_enc_arr_example(mode)
% Plot examples of sh_enc_arr.m and array designs

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%mode:          String, array type {'uni_circ', 'spis_circ', 'spis_custom'}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: 

%plot_sh_enc_arr_example('uni_circ');     % N = 3 circular phase-array
%plot_sh_enc_arr_example('spis_circ');    % N = 3 large piston circular array near equator of spherical baffle  
%plot_sh_enc_arr_example('spis_custom');  % N = 3 micro piston circular array near pole of spherical baffle  

arguments
    mode (1,:) char {mustBeMember(mode, {'uni_circ', 'spis_circ', 'spis_custom'})} = 'uni_circ';
end

if strcmp(mode, 'uni_circ')

    % Sample usage: Design MVDR beamformer over circular omni-directivity transducer array
    % Note: Uniform line and circular arrays yield for phased arrays have real responses for all and even N respectively
    
    Fs = 16000;
    M = 256;
    N = 3;
    radius = 0.1;
    delay_dist = 1;
    max_odr = 36;
    mode = 'uni_circ';
    freq = logspace(log10(20), log10(Fs), M + 1); freq = freq(1:end-1);
    
    theta = pi/2;
    phi = 0;
    kappa = inf;
    mode_mvdr = 'uniform';
    
    C_arr = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, 'enable_disp', true, 'disp_dB_lim', [-3, 3]);
    bf_mvdr_sh(C_arr, false, freq, theta, phi, kappa, mode_mvdr, 'enable_disp', true, 'disp_dB_lim', [-48, 14]);
    
    [C_arr_noncenter, center] = sh_enc_arr(0, mode, N, radius, delay_dist, freq, 'noncentered', true, 'enable_disp', true, 'disp_dB_lim', [-3, 3]);
    bf_mvdr_sh(C_arr_noncenter, false, freq, theta, phi, kappa, mode_mvdr, 'center', center, 'wavnum', 2 * pi * freq / 343, 'enable_disp', true, 'disp_dB_lim', [-48, 14]);

elseif strcmp(mode, 'spis_circ')

    Fs = 24000;
    M = 512;
    N = 3;
    radius = 0.1;
    delay_dist = 1;
    freq = logspace(log10(20), log10(Fs), M + 1); freq = freq(1:end-1);
    spis_arc_len = 0.0254;
    
    max_odr = 50;
    [C_arr, center, h_arr_figs] = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, ...
        'spis_arc_len', spis_arc_len, 'enable_disp', true, 'disp_dB_lim', [-36, 12]);

elseif strcmp(mode, 'spis_custom')

    Fs = 24000;
    M = 512;
    N = 3;
    radius = 0.1;
    delay_dist = 1;
    freq = logspace(log10(20), log10(Fs), M + 1); freq = freq(1:end-1);
    
    spis_arc_len = 0.004;
    
    custom_theta_phi = [deg2rad(10) * ones(1,N); deg2rad(0:(360/N):(360 - 360/N))];

    max_odr = 50;
    [C_arr, center, h_arr_figs] = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, ...
        'spis_arc_len', spis_arc_len, 'custom_theta_phi', custom_theta_phi, ...
        'enable_disp', true, 'disp_dB_lim', [-36, 12]);

    % Design MVDR beam with R, uniform density, on-axis steering angle 
    theta_steer = deg2rad(90);
    phi_steer = deg2rad(0);
    kappa = inf;
    mode_mvdr = 'uniform';
    dB_lim = [-40, 20];
    [~, ~, ~, ~, h_figs_MVDR_uniform_onaxis] = bf_mvdr_sh(C_arr, false, freq, theta_steer, phi_steer, kappa, mode_mvdr, ...
        'enable_disp', true, 'disp_dB_lim', dB_lim);

    % Design MVDR beam with R, uniform density, steer between elements
    theta_steer = deg2rad(90);
    phi_steer = deg2rad(60);
    kappa = inf;
    mode_mvdr = 'uniform';
    [~, ~, ~, ~, h_figs_MVDR_uniform_60deg] = bf_mvdr_sh(C_arr, false, freq, theta_steer, phi_steer, kappa, mode_mvdr, ...
        'enable_disp', true, 'disp_dB_lim', dB_lim);

    % Design MVDR beam with R, custom density nulls broad 90 degree azimuth region 
    theta_steer = deg2rad(90);
    phi_steer = deg2rad(0);
    kappa = inf;
    mode_mvdr = 'density';
    R_pdf = sh_pdf_preset('East', 'SqExp', 6, false, 'enable_disp', true);

    [~, ~, ~, ~, h_figs_MVDR_density_east_broad] = bf_mvdr_sh(C_arr, false, freq, theta_steer, phi_steer, kappa, mode_mvdr, ...
        'R_pdf', R_pdf, ...
        'enable_disp', true, 'disp_dB_lim', dB_lim);

    % Design MVDR beam with R, custom density nulls narrow 90 degree azimuth region 
    theta_steer = deg2rad(90);
    phi_steer = deg2rad(0);
    kappa = inf;
    mode_mvdr = 'density';
    R_pdf = sh_nrm(sh_enc_rbf('SqExp', 30, pi/2, deg2rad(90), 0.15, false), 'Sum');

    [~, ~, ~, ~, h_figs_MVDR_density_east_narrow] = bf_mvdr_sh(C_arr, false, freq, theta_steer, phi_steer, kappa, mode_mvdr, ...
        'R_pdf', R_pdf, ...
        'enable_disp', true, 'disp_dB_lim', dB_lim);

    % Design MVDR beam with R, custom density nulls 10 degree azimuth region
    theta_steer = deg2rad(90);
    phi_steer = deg2rad(0);
    kappa = inf;
    mode_mvdr = 'density';
    R_pdf = sh_nrm(sh_enc_rbf('SqExp', 30, pi/2, deg2rad(10), 0.15, false), 'Sum');

    [~, ~, ~, ~, h_figs_MVDR_density_contra] = bf_mvdr_sh(C_arr, false, freq, theta_steer, phi_steer, kappa, mode_mvdr, ...
        'R_pdf', R_pdf, ...
        'enable_disp', true, 'disp_dB_lim', dB_lim);

    % Design MVDR beam with R, nulls 90 degree azimuth direction, max cond = 1e2
    theta_steer = deg2rad(90);
    phi_steer = deg2rad(0);
    kappa = 1e2;
    mode_mvdr = 'direction';
    R_theta_phi = deg2rad([90, 90]);

    [~, ~, ~, ~, h_figs_MVDR_dir_kappa_1e2] = bf_mvdr_sh(C_arr, false, freq, theta_steer, phi_steer, kappa, mode_mvdr, ...
        'R_theta_phi', R_theta_phi, ...
        'enable_disp', true, 'disp_dB_lim', dB_lim);

    % Design MVDR beam with R, nulls 90 degree azimuth direction, max cond = 1e8
    theta_steer = deg2rad(90);
    phi_steer = deg2rad(0);
    kappa = 1e8;
    mode_mvdr = 'direction';
    R_theta_phi = deg2rad([90, 90]);

    [~, ~, ~, ~, h_figs_MVDR_dir_kappa_1e8] = bf_mvdr_sh(C_arr, false, freq, theta_steer, phi_steer, kappa, mode_mvdr, ...
        'R_theta_phi', R_theta_phi, ...
        'enable_disp', true, 'disp_dB_lim', dB_lim);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Export figures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out_dir = 'figs/figs_bf/arr';
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    exportgraphics(h_figs_MVDR_uniform_onaxis{1}, fullfile(out_dir, 'fit_pist_MVDR_N_3_onaxis_bp.png' ));
    exportgraphics(h_figs_MVDR_uniform_onaxis{2}, fullfile(out_dir, 'fit_pist_MVDR_N_3_onaxis_bs.png' ));

    exportgraphics(h_figs_MVDR_uniform_60deg{1}, fullfile(out_dir, 'fit_pist_MVDR_N_3_60deg_bp.png' ));
    exportgraphics(h_figs_MVDR_uniform_60deg{2}, fullfile(out_dir, 'fit_pist_MVDR_N_3_60deg_bs.png' ));

    exportgraphics(h_figs_MVDR_density_east_broad{1}, fullfile(out_dir, 'fit_pist_MVDR_N_3_estbrd_bp.png' ));
    exportgraphics(h_figs_MVDR_density_east_broad{2}, fullfile(out_dir, 'fit_pist_MVDR_N_3_estbrd_bs.png' ));

    exportgraphics(h_figs_MVDR_density_east_narrow{1}, fullfile(out_dir, 'fit_pist_MVDR_N_3_estnrw_bp.png' ));
    exportgraphics(h_figs_MVDR_density_east_narrow{2}, fullfile(out_dir, 'fit_pist_MVDR_N_3_estnrw_bs.png' ));

    exportgraphics(h_figs_MVDR_density_contra{1}, fullfile(out_dir, 'fit_pist_MVDR_N_3_contra_bp.png' ));
    exportgraphics(h_figs_MVDR_density_contra{2}, fullfile(out_dir, 'fit_pist_MVDR_N_3_contra_bs.png' ));

    exportgraphics(h_figs_MVDR_dir_kappa_1e2{1}, fullfile(out_dir, 'fit_pist_MVDR_N_3_arrg1e2_bp.png' ));
    exportgraphics(h_figs_MVDR_dir_kappa_1e2{2}, fullfile(out_dir, 'fit_pist_MVDR_N_3_arrg1e2_bs.png' ));

    exportgraphics(h_figs_MVDR_dir_kappa_1e8{1}, fullfile(out_dir, 'fit_pist_MVDR_N_3_arrg1e8_bp.png' ));
    exportgraphics(h_figs_MVDR_dir_kappa_1e8{2}, fullfile(out_dir, 'fit_pist_MVDR_N_3_arrg1e8_bs.png' ));

    ;

else
    error('Unknown mode');
end