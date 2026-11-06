function h_fig = bf_plt_sh(C_arr, is_real, W, freq, mode, options)
% Plot spherical harmonic beamformer's beampattern response and statistics

% Detailed beampatterns at select frequencies can be specified in options.disp_freq

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C_arr:     [(L_C + 1)^2 X N x M]   SH expansion of array, L_C max-order, N channels, M frequencies
%is_real:   Logical, if true, evaluate real SH

%W:         [N x M]     Beamformer weights (conjugated)
%freq:      [1 x M]     Frequency (Hz)

%mode:      String, display modes {'horizontal', 'median', 'frontal', 'stats'}
%               'horizontal':   Beampattern in z = 0 plane
%               'median':       Beampattern in y = 0 plane
%               'frontal:       Beampattern in x = 0 plane
%               'stats:         Statistics / frequency (white-noise gain, directivity index)

%options:               struct
%options.center:        [] or [N x 3]   Center-of-expansion coordinates xyz per channel 
%options.wavnum:        [] or [1 x M]   Angular wavenumber kappa (e.g. kappa = 2 * pi * hz / c for propagation velocity c)

%options.resolution:    [1 x 2] Number of sampling points for [theta, phi], must be positive integer
%options.FaceColor:     String, plot interpolation mode {'flat', 'interp'}
%options.fontsize:      Scalar, font size positive integer

%options.disp_mag:      Logical, display magnitude response plots if true
%options.disp_phase:    Logical, display phase response plots if true
%options.disp_cb:       Logical, display color plot if true

%options.disp_title:    Logical, display title label if true
%options.disp_xlabel:   Logical, display xaxis label if true
%options.disp_ylabel:   Logical, display yaxis label if true
%options.disp_xticks:   Logical, display xaxis ticks if true
%options.disp_yticks:   Logical, display yaxis ticks if true

%options.disp_freq:     [1 x M_f]   Display beampatterns at select frequencies if non-empty

%options.dB_lim:        [1 x 2] Magnitude dB range [min, max]
%options.deg_lim:       [1 x 2] Phase degree range [min, max]

%options.title_name:     String, custom title name
%options.title_name_override:   String, custom title name override all names

%options.fig_pos        [1 x 2] [x, y] pixels
%options.fig_size       [1 x 2] [width, height] pixels

%options.stats_DI_theta_phi:    [1 x 2] Steering direction (co-latitude, azimuth radians)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h_fig:             Figure handle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample usage:     Line array, unity beamformer weights (broadside beam)

% Fs = 24000;
% M = 512;
% N = 4;
% %radius = 0.0254; % 1 inch
% radius = 0.0254 * 4; % 4 inch
% delay_dist = 1;
% max_odr = 30;
% mode = 'uni_line';
% freq = logspace(log10(20), log10(Fs), M + 1); freq = freq(1:end-1);
% wavnum = 2 * pi * freq / 343;

% C_arr = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, 'enable_disp', true, 'disp_dB_lim', [-3, 3]);
% W = ones(N, M);
% bf_plt_sh(C_arr, false, W, freq, 'horizontal', 'dB_lim', [-48, 14], 'title_name', 'Centered Array');

% % Reduce max-order expansion to 0 for noncenter origin
% [C_arr_noncenter, center] = sh_enc_arr(0, mode, N, radius, delay_dist, freq, 'noncentered', true, 'enable_disp', true, 'disp_dB_lim', [-3, 3]);

% bf_plt_sh(C_arr_noncenter, false, W, freq, 'horizontal', 'center', center, 'wavnum', wavnum, 'dB_lim', [-48, 14], 'title_name', 'Non-centered Array');
% bf_plt_sh(C_arr_noncenter, false, W, freq, 'median', 'center', center, 'wavnum', wavnum, 'dB_lim', [-48, 14], 'title_name', 'Non-centered Array');
% bf_plt_sh(C_arr_noncenter, false, W, freq, 'frontal', 'center', center, 'wavnum', wavnum, 'dB_lim', [-48, 14], 'title_name', 'Non-centered Array');
% bf_plt_sh(C_arr_noncenter, false, W, freq, 'stats', 'center', center, 'wavnum', wavnum, 'dB_lim', [-48, 14], 'title_name', 'Non-centered Array');

arguments
    C_arr (:,:,:) double {coder.mustBeComplex} = complex(0)
    is_real (1,1) logical = false;

    W (:,:) double  {coder.mustBeComplex} = complex(0);
    freq (1,:) double = 0;

    mode (1,:) char {mustBeMember(mode, {'horizontal', 'median', 'frontal', 'stats'})} = 'horizontal';

    options.center (:,3) double = [];
    options.wavnum (1,:) double = [];

    options.resolution (1,2) double {mustBePositive, mustBeInteger} = [400 400];
    options.FaceColor (1,:) char {mustBeMember(options.FaceColor, {'flat', 'interp'}) } = 'flat';
    options.fontsize (1,1) double {mustBePositive, mustBeInteger} = 16;

    options.disp_mag (1,1) logical = true;
    options.disp_phase (1,1) logical = true;
    options.disp_cb (1,1) logical = true;

    options.disp_title  (1,1) logical = true;
    options.disp_xlabel (1,1) logical = true;
    options.disp_ylabel (1,1) logical = true;
    options.disp_xticks (1,1) logical = true;
    options.disp_yticks (1,1) logical = true;
    
    options.disp_freq (1,:) double {mustBeNonnegative} = [];

    options.dB_lim (1,2) double  = [-inf inf];
    options.deg_lim (1,2) double = [-180 180];

    options.title_name (1,:) char = [];
    options.title_name_override (1,:) char = [];

    options.fig_pos  (1,2) double = [100, 100];
    options.fig_size (1,2) double {mustBePositive} = [800, 600];

    options.stats_DI_theta_phi (1,2) double = [pi/2, 0];

end

[N_C, N, M] = size(C_arr);

if isempty(options.center)
    center_wavnum = [];
else
    center_wavnum = zeros([N, M, 4]);
    for n = 1:N
        for m = 1:M
            center_wavnum(n, m, :) = [options.center(n, :), options.wavnum(m) ];
        end
    end
    center_wavnum = reshape(center_wavnum, [N * M, 4]);
end

if any(strcmp(mode, {'horizontal', 'median', 'frontal'}))

    if strcmp(mode, 'horizontal')
    
        N_dir = options.resolution(2);
        phi = linspace(-pi, pi, N_dir)';
        theta = pi/2 * ones(N_dir, 1);
        planar_label = 'Azimuth (Degrees)';
        dir_deg = rad2deg(phi);

    elseif strcmp(mode, 'median')

        N_dir = options.resolution(1);
        phi = zeros(N_dir, 1);
        theta = pi/2 + linspace(-pi, pi, N_dir)';
        planar_label = 'Elevation (Degrees)';
        dir_deg = rad2deg(pi/2 - theta);

    elseif strcmp(mode, 'frontal')

        N_dir = options.resolution(1);
        phi = pi/2 * ones(N_dir, 1);
        theta = pi/2 + linspace(-pi, pi, N_dir)';
        % First dir is along -y, then increases in elevation
        planar_label = 'Elevation (Degrees)';
        dir_deg = rad2deg(pi/2 - theta);

    else

        error('Unknown mode');

    end

    % Compute array responses
    f_arr = reshape(sh_dec(reshape(C_arr, [N_C, N * M]), theta, phi, is_real, 'eval', center_wavnum), [N_dir, N, M] ); % [N_dir x N x M]

    % Compute beampattern responses
    f = complex(zeros([N_dir, M]));
    for m = 1:M
        f(:, m) = f_arr(:, :, m) * conj(W(:, m)); % [N_dir x 1]        
    end

    h_fig = figure;
    h_fig.Position = [options.fig_pos, options.fig_size];
    if options.disp_mag && options.disp_phase
        h_tile = tiledlayout(2,1);
    elseif options.disp_mag || options.disp_phase
        h_tile = tiledlayout(1,1);
    else
        error('options.disp_mag or options.disp_phase must be true');
    end

    if options.disp_mag % Plot magnitude response
        nexttile;

        assert(numel(freq) == M, 'freq size mismatch C');

        if isempty(options.disp_freq) % Display all frequencies
        
            h_pc = pcolor(freq, dir_deg, mag2db(abs(f)) );
            xscale("log");
    
           % Plot grid
            freq_lines = [10:10:90, 100:100:900, 1000:1000:9000, 10000:10000:20000];
            freq_lines(freq_lines < min(freq)) = [];
            freq_lines(freq_lines > max(freq)) = [];
            for i = 1:numel(freq_lines)
                xline(freq_lines(i), '--');
            end
            yline(0, 'r--'); 
            yline(90, 'm--');
            yline(-90, 'm--');
            yline(45, 'g--');
            yline(-45, 'g--');
    
            set(h_pc, 'EdgeColor', 'none');
            set(h_pc, 'FaceColor', options.FaceColor);
    
            if options.disp_xlabel
                xlabel('Frequency (Hz)', 'fontsize', options.fontsize);            
            end
            if ~options.disp_xticks
                set(gca,'XTick', []);
            end
            if options.disp_ylabel
                ylabel(planar_label, 'fontsize', options.fontsize);
                %ylabel('Azimuth', 'fontsize', options.fontsize);
            end
            if ~options.disp_yticks
                set(gca,'YTick', []);
            end
    
            clim([options.dB_lim(1), options.dB_lim(2)]);
            if options.disp_cb
                h_cb = colorbar;
                ylabel(h_cb, 'dB', 'fontsize', options.fontsize);
            end
            %xscale log;

        else % Display select frequencies

            xline(0, 'r--', 'linewidth', 1.5, 'HandleVisibility', 'off'); hold on;
            xline(90, 'm--', 'linewidth', 1.5, 'HandleVisibility', 'off');
            xline(-90, 'm--', 'linewidth', 1.5, 'HandleVisibility', 'off');
            xline(45, 'g--', 'linewidth', 1.5, 'HandleVisibility', 'off');
            xline(-45, 'g--', 'linewidth', 1.5, 'HandleVisibility', 'off');

            M_f = numel(options.disp_freq);
            legend_str = cell(1, M_f);
            for m = 1:M_f
                
                [~, idx_m] = min(abs(freq - options.disp_freq(m))); idx_m = idx_m(1);
                plot( dir_deg, mag2db(abs(f(:, idx_m))), 'linewidth', 2.5); hold on;
                legend_str{m} = num2str(freq(idx_m));

            end
            h_lg = legend(legend_str, 'location', 'best');
            set(h_lg, 'fontsize', options.fontsize - 1); title(h_lg, 'Freq. (Hz)');
            grid on; axis tight;

            if options.disp_xlabel
                xlabel(planar_label, 'fontsize', options.fontsize);
            end
            if ~options.disp_xticks
                set(gca,'XTick', []);
            end
            if options.disp_ylabel
                ylabel('Magnitude (dB)', 'fontsize', options.fontsize);
            end
            if ~options.disp_yticks
                set(gca,'YTick', []);
            end
            
        end

        title_str = 'Magnitude Response';
        if ~isempty(options.title_name)
            title_str = [options.title_name, ' ', title_str];
        end
        if options.disp_title
            if isempty(options.title_name_override)
                title(title_str, 'fontsize', options.fontsize + 1);  
            else
                title(options.title_name_override, 'fontsize', options.fontsize + 1);
            end  
        end

        set(gca, 'fontsize', options.fontsize - 1);

    end

    if options.disp_phase % Plot phase response
        nexttile;

        assert(numel(freq) == M, 'freq size mismatch C');

        if isempty(options.disp_freq) % Display all frequencies
    
            h_pc = pcolor(freq, dir_deg, rad2deg(angle(f)) ); hold on;
            xscale("log");
    
            % Plot grid
            freq_lines = [10:10:90, 100:100:900, 1000:1000:9000, 10000:10000:20000];
            freq_lines(freq_lines < min(freq)) = [];
            freq_lines(freq_lines > max(freq)) = [];
            for i = 1:numel(freq_lines)
                xline(freq_lines(i), '--');
            end
            yline(0, 'r--'); 
            yline(90, 'm--');
            yline(-90, 'm--');
            yline(45, 'g--');
            yline(-45, 'g--');
    
            set(h_pc, 'EdgeColor', 'none');
            set(h_pc, 'FaceColor', options.FaceColor);
    
            if options.disp_xlabel
                xlabel('Frequency (Hz)', 'fontsize', options.fontsize);            
            end
            if ~options.disp_xticks
                set(gca,'XTick', []);
            end
            if options.disp_ylabel
                ylabel(planar_label, 'fontsize', options.fontsize);
            end
            if ~options.disp_yticks
                set(gca,'YTick', []);
            end          
    
            clim([options.deg_lim(1), options.deg_lim(2)]);
            if options.disp_cb
                h_cb = colorbar;
                ylabel(h_cb, 'Degrees', 'fontsize', options.fontsize);
            end
            colormap(gca, hsv);
        
        else % Display select frequencies

            xline(0, 'r--', 'linewidth', 1.5, 'HandleVisibility', 'off'); hold on;
            xline(90, 'm--', 'linewidth', 1.5, 'HandleVisibility', 'off');
            xline(-90, 'm--', 'linewidth', 1.5, 'HandleVisibility', 'off');
            xline(45, 'g--', 'linewidth', 1.5, 'HandleVisibility', 'off');
            xline(-45, 'g--', 'linewidth', 1.5, 'HandleVisibility', 'off');

            M_f = numel(options.disp_freq);
            legend_str = cell(1, M_f);
            for m = 1:M_f
                
                [~, idx_m] = min(abs(freq - options.disp_freq(m))); idx_m = idx_m(1);
                plot( dir_deg, rad2deg(angle(f(:, idx_m))), 'linewidth', 2.5); hold on;
                legend_str{m} = num2str(freq(idx_m));

            end
            h_lg = legend(legend_str, 'location', 'best');
            set(h_lg, 'fontsize', options.fontsize - 1); title(h_lg, 'Freq. (Hz)');
            grid on; axis tight;

            if options.disp_xlabel
                xlabel(planar_label, 'fontsize', options.fontsize);
            end
            if ~options.disp_xticks
                set(gca,'XTick', []);
            end
            if options.disp_ylabel
                ylabel('Phase (Degrees)', 'fontsize', options.fontsize);
            end
            if ~options.disp_yticks
                set(gca,'YTick', []);
            end

        end

        title_str = 'Phase Response';
        if ~isempty(options.title_name)
            title_str = [options.title_name, ' ', title_str];
        end
        if options.disp_title            
            if isempty(options.title_name_override)
                title(title_str, 'fontsize', options.fontsize + 1);  
            else
                title(options.title_name_override, 'fontsize', options.fontsize + 1);
            end
        end
        
        set(gca, 'fontsize', options.fontsize - 1);

    end

    % Set global title
    title(h_tile, ['Beam Pattern Responses (', mode, ' plane)'], 'fontsize', options.fontsize + 2);

elseif strcmp(mode, 'stats')

    % Compute electrical power
    EP_dB = pow2db(sum(W .* conj(W), 1)); % Total electrical power
    CP_dB = pow2db(W .* conj(W));         % Channel power [N x M]

    % Compute directivity index
    R = sh_cov(C_arr, sh_pdf_preset('FwdCenter', 'Uniform', 0), is_real, 'SHT_center', options.center, 'SHT_wavnum', options.wavnum, ...
        'SHT_oversample', 1, 'SHT_maxsample', 1024); % [N x N x M] Spatial covariance integrated over uniform density

    theta = options.stats_DI_theta_phi(1);
    phi = options.stats_DI_theta_phi(2);    
    V = reshape(sh_dec(reshape(C_arr, [N_C, N * M]), theta, phi, is_real, 'eval', center_wavnum), [N, M] ); % [N x M] steering vectors
    DI_dB = zeros([1, M]);
    CDI_dB = zeros([N, M]);
    for m = 1:M
        DI_dB(m) = pow2db( abs(W(:,m)' * V(:,m) )^2 / real( W(:,m)' * R(:,:,m) * W(:,m) ) );
        for n = 1:N
            CDI_dB(n, m) = pow2db( abs(W(n,m)' * V(n,m) )^2 / real( W(n,m)' * R(n,n,m) * W(n,m) ) );
        end
    end

    % Display
    h_fig = figure;
    h_fig.Position = [options.fig_pos, options.fig_size];
    h_tile = tiledlayout(2,1);
        
    % Plot directivity index
    nexttile;
    semilogx(freq, DI_dB, 'b-', 'linewidth', 3); 
    legend_str = {'Array'}; hold on;
    for n = 1:N
        semilogx(freq, CDI_dB(n,:), '-', 'linewidth', 1.5); 
        legend_str{end+1} = ['Ch. ', num2str(n)];
    end
    grid on; axis tight;
    xlabel('Frequency (Hz)', 'fontsize', options.fontsize);
    ylabel('Magnitude (dB)', 'fontsize', options.fontsize);

    h_lg = legend(legend_str, 'location', 'best', 'NumColumns', 6, 'Orientation', 'horizontal');
    set(h_lg, 'fontsize', options.fontsize - 1);
    set(gca, 'fontsize', options.fontsize - 1);
    title('Directivity Index', 'fontsize', options.fontsize + 1);

    % Plot electrical power
    nexttile;
    semilogx(freq, EP_dB, 'r-', 'linewidth', 3); 
    legend_str = {'Total'}; hold on;
    for n = 1:N
        semilogx(freq, CP_dB(n,:), '-', 'linewidth', 1.5); 
        legend_str{end+1} = ['Ch. ', num2str(n)];
    end
    grid on; axis tight;
    xlabel('Frequency (Hz)', 'fontsize', options.fontsize);
    ylabel('Magnitude (dB)', 'fontsize', options.fontsize);
        
    h_lg = legend(legend_str, 'location', 'best', 'NumColumns', 6,  'Orientation', 'horizontal');
    set(h_lg, 'fontsize', options.fontsize - 1);
    set(gca, 'fontsize', options.fontsize - 1);
    title('Electrical Power', 'fontsize', options.fontsize + 1);



    % Global title
    title_str = 'Beamformer Statistics: ';
    if ~isempty(options.title_name)
        title_str = [title_str, ' ', options.title_name];
    end
    if isempty(options.title_name_override)
        title(h_tile, title_str, 'fontsize', options.fontsize + 1);  
    else
        title(h_tile, options.title_name_override, 'fontsize', options.fontsize + 2);
    end
            

else
    error('Unknown mode');
end

