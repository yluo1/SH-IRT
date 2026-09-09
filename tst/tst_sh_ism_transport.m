function tst_sh_ism_transport(mode)
%Test spherical harmonic image-source model (sh_ism.m)
%with PDF transport of echo densities 

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%mode:      String, test case {'barycenter', 'path_dist'}

%options.use_mex:   Logical, if true, use mex sub-functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:

%tst_sh_ism_transport('barycenter');
%tst_sh_ism_transport('path_dist');

arguments
    mode (1,:) char {mustBeMember(mode, {'barycenter', 'path_dist'})} = 'barycenter';
end

max_src_odr = 1;

%max_rec_odr = 4;
max_rec_odr = 8;

is_real = true;
%is_real = false;

T = 0.155;                    %Time (Sec) duration

%l   = [5    6       4];     %Room size
l   = [7    8       5];     %Room size

r   = [0.1  0.3    0.5];    %Receiver

    
% gamma_pos = [0.9, 0.9, 0.9];
% gamma_neg = [0.9, 0.9, 0.9];

% gamma_pos = ones(1,3) * 0.8;
% gamma_neg = ones(1,3) * 0.8;

% gamma_pos = ones(1,3) * 0.5;
% gamma_neg = ones(1,3) * 0.5;

gamma_pos = ones(1,3) * (1/3);
gamma_neg = ones(1,3) * (1/3);

% gamma_pos = ones(1,3) * 0.25;
% gamma_neg = ones(1,3) * 0.25;

% gamma_pos = [0.9 0.8 0.85] * 0.25;
% gamma_neg = [0.92 0.85 0.9] * 0.25;

% gamma_pos = ones(1,3) * 0.1;
% gamma_neg = ones(1,3) * 0.1;


%s_C = [0.2  -2.9    1.9];   %Source location, left of receiver
%s_C = [2 -2.5    1.9];    %Source location, upper left corner
%s_C = [2    -2.5   0];      %Source location,  left corner
s_C = [1    -1   0];      %Source location,  left forward
%s_C = [1    -1   1];      %Source location,  left forward elevated

% s_D = [2  2.5   -1.9];   %Source location, bottom right corner
%s_D = [2     2.5   0];      %Source location,  right corner
s_D = [1     1   0];      %Source location,  right forward
%s_D = [1     1   -1];      %Source location,  right forward lowered

if strcmp(mode, 'barycenter')

    s_E = (s_C + s_D) / 2;       %Source average location
    
    % Generate RIR density
    varargin_ism = {'mode_enc', 'proj_msq_mex', 'jitter_coord_bnd', [1 1] * 1e-3};
    [C, t_C] = sh_ism(max_src_odr, max_rec_odr, is_real, T, s_C, r, l, gamma_pos, gamma_neg, varargin_ism{:});
    [D, t_D] = sh_ism(max_src_odr, max_rec_odr, is_real, T, s_D, r, l, gamma_pos, gamma_neg, varargin_ism{:});
    [E, t_E] = sh_ism(max_src_odr, max_rec_odr, is_real, T, s_E, r, l, gamma_pos, gamma_neg, varargin_ism{:});
    
    %plot_RIR(squeeze(C(1, 1, :)), plot_RIR_opts('clim', [-120, -60], 'win_size', 512, 'spec_scale', 'linear'));
        
    %Integrated density in the specified time interval
    t_lim = [0, inf];      % All time
    %t_lim = [0, 0.05];      % Direct + Early reflections
    %t_lim = [0.05, inf];   % Late reflections
    
    t_C_idx = t_C >= min(t_lim) & t_C <= max(t_lim);
    t_D_idx = t_D >= min(t_lim) & t_D <= max(t_lim);
    t_E_idx = t_E >= min(t_lim) & t_E <= max(t_lim);
    
    C_pdf = sum(squeeze(C(1, :, t_C_idx)), 2);
    C_pdf = sh_nrm(C_pdf, 'Sum');
    
    D_pdf = sum(squeeze(D(1, :, t_D_idx)), 2);
    D_pdf = sh_nrm(D_pdf, 'Sum');
    
    E_pdf = sum(squeeze(E(1, :, t_E_idx)), 2);
    E_pdf = sh_nrm(E_pdf, 'Sum');
    
    [E_pdf_linear] = sh_pdf_transport(C_pdf, D_pdf, 0.5, 'LinearInterp', is_real);
    [E_pdf_geometric] = sh_pdf_transport(C_pdf, D_pdf, 0.5, 'GeometricInterp', is_real);
    [E_pdf_wasserstein] = sh_pdf_transport(C_pdf, D_pdf, 0.5, 'SlicedWass', is_real);
    

    %Plot Horizontal map
    fontsize = 18;
    varargin_sh_plt_horizontal = {'disp_phase', false, 'dB_lim', [-120, -20], 't', t_C, 'disp_xaxis_ker_size', 1024, 'fig_size', [500, 300] * (3/4), 'fontsize', fontsize};
    h_RIR_left      = sh_plt(squeeze(C(1, :, :)), 'horizontal', is_real, varargin_sh_plt_horizontal{:}, 'title_name_override', 'SH-ISM Expansion (t = 0)', 'disp_cb', false);
    h_RIR_center    = sh_plt(squeeze(E(1, :, :)), 'horizontal', is_real, varargin_sh_plt_horizontal{:}, 'title_name_override', 'SH-ISM Expansion (t = 0.5)', 'disp_ylabel', false, 'disp_yticks', false, 'disp_cb', false);
    h_RIR_right     = sh_plt(squeeze(D(1, :, :)), 'horizontal', is_real, varargin_sh_plt_horizontal{:}, 'title_name_override', 'SH-ISM Expansion (t = 1)', 'disp_ylabel', false, 'disp_yticks', false,  'disp_cb', false);

    %Plot density
    varargin_sh_plt_density = {'disp_phase', false, 'dB_lim', [-38, 6], ...
        'fig_pos', [800, 100], 'fig_size', [500, 300] * (3/4), 'fontsize', fontsize};
    
    h_density_left =    sh_plt(C_pdf, 'mercator', is_real, varargin_sh_plt_density{:}, 'title_name_override', 'RIR Start Density Mean', 'disp_cb', false);
    h_density_center =  sh_plt(E_pdf, 'mercator', is_real, varargin_sh_plt_density{:}, 'title_name_override', 'RIR Center Density Mean', 'disp_cb', false, 'disp_ylabel', false, 'disp_yticks', false);
    h_density_right =   sh_plt(D_pdf, 'mercator', is_real, varargin_sh_plt_density{:}, 'title_name_override',  'RIR End Density Mean', 'disp_ylabel', false, 'disp_yticks', false,  'disp_cb', false);
    
    h_arithemtic    = sh_plt(E_pdf_linear, 'mercator', is_real, varargin_sh_plt_density{:}, 'title_name_override', 'Linear Density (t = 0.5)', 'disp_cb', false);
    h_geometric     = sh_plt(E_pdf_geometric, 'mercator', is_real, varargin_sh_plt_density{:}, 'title_name_override', 'Product Density (t = 0.5)', 'disp_cb', false, 'disp_ylabel', false, 'disp_yticks', false);
    h_barycentric   = sh_plt(E_pdf_wasserstein, 'mercator', is_real, varargin_sh_plt_density{:}, 'title_name_override', 'Wasserstein Density (t = 0.5)', 'disp_ylabel', false, 'disp_yticks', false,  'disp_cb', false);
    
    
    %export
    out_dir = 'figs_sh';
    mkdir(out_dir);
    exportgraphics(h_RIR_left{1} , fullfile(out_dir, ['RIR_left', '.png']));
    exportgraphics(h_RIR_right{1} , fullfile(out_dir, ['RIR_right', '.png']));
    exportgraphics(h_RIR_center{1} , fullfile(out_dir, ['RIR_center', '.png']));
    
    exportgraphics(h_density_left{1} , fullfile(out_dir, ['density_left', '.png']));
    exportgraphics(h_density_right{1} , fullfile(out_dir, ['density_right', '.png']));
    exportgraphics(h_density_center{1} , fullfile(out_dir, ['density_center', '.png']));
    
    exportgraphics(h_arithemtic{1} , fullfile(out_dir, ['arithemtic_mean', '.png']));
    exportgraphics(h_geometric{1} , fullfile(out_dir, ['geometric_mean', '.png']));
    exportgraphics(h_barycentric{1} , fullfile(out_dir, ['barycentric_mean', '.png']));

    ;

elseif strcmp(mode, 'path_dist')

    N = 11; %Number of points along line-path
    %N = 3;

    t = linspace(0, 1, N)';

    % Generate RIR density
    varargin_ism = {'mode_enc', 'proj_msq_mex', 'jitter_coord_bnd', [1 1] * 1e-3};
    [C, t_C] = sh_ism(max_src_odr, max_rec_odr, is_real, T, s_C, r, l, gamma_pos, gamma_neg, varargin_ism{:});
    [D, t_D] = sh_ism(max_src_odr, max_rec_odr, is_real, T, s_D, r, l, gamma_pos, gamma_neg, varargin_ism{:});

    %Integrated density in the specified time interval
    t_lim = [0, inf];      % All time
    %t_lim = [0, 0.05];      % Direct + Early reflections
    %t_lim = [0.05, inf];   % Late reflections
    
    t_C_idx = t_C >= min(t_lim) & t_C <= max(t_lim);
    t_D_idx = t_D >= min(t_lim) & t_D <= max(t_lim);
    
    s_E = (1-t) * s_C + t * s_D; %Source position along path

    C_pdf = sum(squeeze(C(1, :, t_C_idx)), 2);
    C_pdf = sh_nrm(C_pdf, 'Sum');
    
    D_pdf = sum(squeeze(D(1, :, t_D_idx)), 2);
    D_pdf = sh_nrm(D_pdf, 'Sum');

    %Track distances
    d_linear = zeros(N, 1);
    d_geometric = zeros(N, 1);
    d_wasserstein = zeros(N, 1);

    %Iterate over path
    for n = 1:N

        % Generate RIR density along path
        [E, t_E] = sh_ism(max_src_odr, max_rec_odr, is_real, T, s_E(n,:), r, l, gamma_pos, gamma_neg, varargin_ism{:});
        t_E_idx = t_E >= min(t_lim) & t_E <= max(t_lim);

        E_pdf = sum(squeeze(E(1, :, t_E_idx)), 2);
        E_pdf = sh_nrm(E_pdf, 'Sum');

        [E_pdf_linear]      = sh_pdf_transport(C_pdf, D_pdf, t(n), 'LinearInterp', is_real);
        [E_pdf_geometric]   = sh_pdf_transport(C_pdf, D_pdf, t(n), 'GeometricInterp', is_real);
        [E_pdf_wasserstein] = sh_pdf_transport(C_pdf, D_pdf, t(n), 'SlicedWass', is_real);

        d_linear(n)         = sh_pdf_dist(E_pdf, E_pdf_linear, 'SlicedWasserstein2', is_real);
        d_geometric(n)      = sh_pdf_dist(E_pdf, E_pdf_geometric, 'SlicedWasserstein2', is_real);
        d_wasserstein(n)    = sh_pdf_dist(E_pdf, E_pdf_wasserstein, 'SlicedWasserstein2', is_real);
    end

    ;

    %Plotting
    fontsize = 18;
    h_fig = figure;
    h_fig.Position = [100, 100, [550, 300] * (3/4)];
    plot(t, d_linear, 'o-', t, d_geometric, 's-', t, d_wasserstein, '*-', 'linewidth', 1.5, 'MarkerSize', 12);
    grid on; axis tight;
    xlabel('Normalized Path Time $t$', 'fontsize', fontsize, 'interpreter', 'latex');
    ylabel('$SSW_2^2(\bf{D}_t, \bf{D}_{t*})$', 'fontsize', fontsize, 'interpreter', 'latex');
    ylim([-inf, max([max(d_linear), max(d_geometric), max(d_wasserstein)]) * 1.1]);
    title('Wasserstein Distances across RIRs', 'fontsize', fontsize + 1);
    set(gca, 'fontsize', fontsize - 1);
    h_lg = legend('Linear', 'Product', 'Wasser.', 'location' ,'northwest');
    set(h_lg, 'fontsize', fontsize - 1);

    %Plot path
    fontsize = 20;
    h_fig_path = figure; 
    h_fig_path.Position = [100, 100, [350, 225] * 1];
    rectangle('Position', [-l(1)/2, -l(2)/2, l(1), l(2) ], 'Linewidth', 1.5 ); hold on;
    plot(s_E(:, 1), s_E(:, 2), 'bx-', r(1), r(2), '*', 'linewidth', 1.5, 'MarkerSize', 12); hold on;
    text(s_E(1, 1) + 0.15, s_E(1, 2) - 2.25, 't = 0 \rightarrow', 'HorizontalAlignment','left', 'fontsize', fontsize - 3)
   % text(s_E(floor(N/2), 1) + 1, s_E(floor(N/2), 2) + 0.35, 't = 0.5', 'HorizontalAlignment', 'center', 'fontsize', fontsize - 3)
    text(s_E(end, 1) + 0.15, s_E(end, 2) + 2.25, '\leftarrow t = 1', 'HorizontalAlignment','right', 'fontsize', fontsize - 3)
    view(-90, 90)
    set(gca, 'YDir','reverse')
    grid on; 
    xlabel('X Dim. (Meters)', 'fontsize', fontsize);
    ylabel('Y Dim. (Meters)', 'fontsize', fontsize);
    title('Room Config.', 'fontsize', fontsize);
    h_lg_path = legend('Source Path', 'Receiver');
    h_lg_path.Position = [0.325, 0.275, 0.4, 0.2];
    set(h_lg_path, 'fontsize', fontsize - 1);
    set(gca, 'fontsize', fontsize - 1);


    %Export
    out_dir = 'figs_sh';
    mkdir(out_dir);
    exportgraphics(h_fig , fullfile(out_dir, ['RIR_Wasserstein_dists', '.png']));
    exportgraphics(h_fig_path , fullfile(out_dir, ['RIR_Wasserstein_room', '.png']));
    ;

else
    error('Unknown mode');
end