function h_fig = plot_sh_irt_logo(mode, options)
% Generate SH-IRT logos

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%mode:      String, logo name {'polar_hplane'}

%options:                   struct
%options.enable_export:     Logical, if true, export figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h_fig:     Handle to figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample usage: Generate sample logos

% plot_sh_irt_logo('polar_hplane')

arguments
    mode (1,:) char {mustBeMember(mode, {'polar_hplane'})} = 'polar_hplane';
    options.enable_export (1,1) logical = false;
end

if strcmp(mode, 'polar_hplane')
    
    N = 1024; % Number of angles
    theta = pi/2 * ones(N, 1);
    phi = linspace(0, 2*pi, N)';
    is_real = false;

    %max_odr = 8;
    max_odr = 9;
    %max_odr = 10;

    Y = sh_val(max_odr, theta, phi, is_real); %[N x (max_odr + 1)^2]
    
    % Plot
    h_fig = figure;
    h_fig.Position = [100, 100, 600, 600];

    linewidth = 3;
    %linewidth = 4;
    
    alpha = 0.075;

    colororder("sail") % 1
%    colororder("earth") % 3
%    colororder("reef") % 2
 %   colororder("gem12") % 4
    
    for l = 0:max_odr
        for m = -l:l
            idx = ((l^2) + 1) + m + l;
            ax = gca;
            set(ax, 'ColorOrderIndex', mod(idx-1, size(gca().ColorOrder, 1))+1 );
            color_lm = [ax.ColorOrder(ax.ColorOrderIndex, :), alpha]; % Color with alpha transparency

            h_p = polarplot(phi, real(Y(:, idx)), '-', 'linewidth',  linewidth);  hold on;
            h_p.Color = color_lm';
            h_p = polarplot(phi, imag(Y(:, idx)), '.-', 'linewidth', linewidth);  hold on;
            h_p.Color = color_lm';
            ;
        end
    end

    % Remove labels and ticks
    ax = gca;
    ax.RTickLabel = [];
    ax.ThetaTickLabel = [];

    % Title
    text(pi/2, 0.65, 0, {'Spherical Harmonic'}, 'fontsize', 36, ...
        'HorizontalAlignment', 'center', 'FontName', 'Andale Mono');
    text(-pi/2, 0.6, 0, {'Impulse Response Tools'}, 'fontsize', 36, ...
        'HorizontalAlignment', 'center', 'FontName', 'Andale Mono');
    
    annotation('rectangle',[0 0 1 1],'Color','w');
    
else
    error('Unknown mode');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.enable_export && coder.target("MATLAB")

    out_dir = 'figs/logo';
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    exportgraphics(h_fig, fullfile(out_dir, ['logo_', mode, '.png'] ), 'Resolution', 300);
end