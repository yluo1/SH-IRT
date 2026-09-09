function tst_sh_ism(mode)
%Test spherical harmonic image-source model (sh_ism.m)

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%mode:      String, test case {'direct_only', 'single_tap', 'multi_tap'}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:

%tst_sh_ism('single_tap');
%tst_sh_ism('multi_tap');

arguments
    mode (1,:) char {mustBeMember(mode, {'direct_only', 'single_tap', 'multi_tap'})} = 'multi_tap';
end

if strcmp(mode, 'direct_only')

    max_src_odr = 1;
    max_rec_odr = 1;
    is_real = true;
    %is_real = false;
    
    T = 0.2;
    
    s = [2 0 1];
    r = [1 0 1];
    l = [5 6 3];
    

    gamma_pos = [0, 0, 0];
    gamma_neg = [0, 0, 0];
    
   
    [C, t] = sh_ism(max_src_odr, max_rec_odr, is_real, T, s, r, l, gamma_pos, gamma_neg, 'jitter_coord_bnd', [1 1] * 1e-3);
    
    plot_RIR(squeeze(C(1, 1, :)), plot_RIR_opts('clim', [-120, -60], 'win_size', 512, 'spec_scale', 'linear'));
    
    % plot_RIR(squeeze(C(2, 1, :)));
    % plot_RIR(squeeze(C(3, 1, :)));
    %plot_RIR(squeeze(C(4, 1, :)));
    %plot_RIR(squeeze(C(4, 4, :)));
    
    sh_plt(squeeze(C(1, :, :)), 'horizontal', is_real, 'disp_phase', false, 'dB_lim', [-120, -10], 't', t, 'disp_xaxis_ker_size', 1024);
    sh_plt(squeeze(C(1, :, ceil( 0.0029 * 48000)))', 'mercator', is_real);

    
elseif    strcmp(mode, 'single_tap')

    max_src_odr = 1;
    max_rec_odr = 5;
    is_real = true;
    %is_real = false;
    
    T = 0.2;
    
    s = [2 0 1];
    r = [1 0 1];
    l = [5 6 3];
    

    gamma_pos = [0.8, 0.7, 0.5];
    gamma_neg = [0.9, 0.6, 0.5];
    
   
    [C, t] = sh_ism(max_src_odr, max_rec_odr, is_real, T, s, r, l, gamma_pos, gamma_neg, 'jitter_coord_bnd', [1 1] * 1e-3);
    
    plot_RIR(squeeze(C(1, 1, :)), plot_RIR_opts('clim', [-120, -60], 'win_size', 512, 'spec_scale', 'linear'));
    

    sh_plt(squeeze(C(1, :, :)), 'horizontal', is_real, 'disp_phase', false, 'dB_lim', [-120, -10], 't', t, 'disp_xaxis_ker_size', 1024);
    sh_plt(squeeze(C(1, :, ceil( 0.0029 * 48000)))', 'mercator', is_real);

    
elseif strcmp(mode, 'multi_tap')

    max_src_odr = 1;
    max_rec_odr = 5;
    is_real = true;
    %is_real = false;
    
    T = 0.2;
    
    s = [2 0 1];
    r = [1 0 1];
    l = [5 6 3];
    
    % gamma_pos = [0.0, 0.0, 0.0];
    % gamma_neg = [0.0, 0.0, 0.0];
    
    % gamma_pos = [0.8, 0.7, 0.5];
    % gamma_neg = [0.9, 0.6, 0.5];
    
    gamma_pos = [0.8, 0.7, 0.5; 0.2, 0.1, 0.3];
    gamma_neg = [0.9, 0.6, 0.5; 0.1, 0.1, 0.2];
    
    
    [C, t] = sh_ism(max_src_odr, max_rec_odr, is_real, T, s, r, l, gamma_pos, gamma_neg, 'jitter_coord_bnd', [1 1] * 1e-3);
    
    plot_RIR(squeeze(C(1, 1, :)), plot_RIR_opts('clim', [-120, -60], 'win_size', 512, 'spec_scale', 'linear'));
    
    % plot_RIR(squeeze(C(2, 1, :)));
    % plot_RIR(squeeze(C(3, 1, :)));
    %plot_RIR(squeeze(C(4, 1, :)));
    %plot_RIR(squeeze(C(4, 4, :)));
    
    sh_plt(squeeze(C(1, :, :)), 'horizontal', is_real, 'disp_phase', false, 'dB_lim', [-120, -10], 't', t, 'disp_xaxis_ker_size', 1024);
    sh_plt(squeeze(C(1, :, ceil( 0.0029 * 48000)))', 'mercator', is_real);

end
