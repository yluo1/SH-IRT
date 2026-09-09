function h = sh_plt(C, mode, is_real, options)
%Plot spherical harmonic (SH) expansions

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:                 [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)

%mode:              String, display mode
%                   'mercator'          
%                   'horizontal'

%is_real:           Logical, if true, evaluate real SH

%options:              struct

%options.hz:            [1 x M]  Frequencies in hz ([] to ignore)
%options.t:             [1 x M]  Time in seconds ([] to ignore)

%options:               struct
%options.resolution:    [1 x 2] Number of sampling points for [theta, phi], must be positive integer
%options.FaceColor:     String, plot interpolation mode {'flat', 'interp'}
%options.fontsize:      Scalar, font size positive integer

%options.rot_theta_zyx  [1 x 3] Rotate expansion by radians along +z, +y, +x (yaw, pitch, roll) axes

%options.disp_mag:      Logical, display magnitude response plots if true
%options.disp_phase:    Logical, display phase response plots if true
%options.disp_cb:       Logical, display color plot if true

%options.disp_title:    Logical, display title label if true
%options.disp_xlabel:   Logical, display xaxis label if true
%options.disp_ylabel:   Logical, display yaxis label if true
%options.disp_xticks:   Logical, display xaxis ticks if true
%options.disp_yticks:   Logical, display yaxis ticks if true

%options.dB_lim:        [1 x 2] Magnitude dB range [min, max]
%options.deg_lim:       [1 x 2] Phase degree range [min, max]

%options.disp_xaxis_ker_size: Size of Hann kernel for increasing solution of xaxis, positive

%options.title_name:     String, custom title name
%options.title_name_override:   String, custom title name override all names

%options.fig_pos        [1 x 2] [x, y] pixels
%options.fig_size       [1 x 2] [width, height] pixels

%options.disp_theta_phi:    [N x 2] (theta, phi). If non-empty [], display coordinates in 'mercator' mode
%options.disp_theta_phi_markersize:   Marker size for displayed theta, phi coordinates 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h:                 [1 x *] Cell array of figure handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Display SH radial basis functions
% P = 3;
% theta = pi/2;
% phi = 0;
% ell = 1;
% is_real = false;
% C = sh_enc_rbf('SqExp', max_odr, theta, phi, ell, is_real, false);

% sh_plt(C, 'mercator', is_real);
% sh_plt(C, 'mercator', is_real, 'rot_theta_zyx', deg2rad([90 0 0]));

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
    mode (1,:) char {mustBeMember(mode, {'mercator', 'horizontal'})} = 'mercator';
    is_real (1,1) logical = false;

    options.hz (1,:) double = [];
    options.t (1,:) double = [];

    options.resolution (1,2) double {mustBePositive, mustBeInteger} = [200 400];
    options.FaceColor (1,:) char {mustBeMember(options.FaceColor, {'flat', 'interp'}) } = 'flat';
    options.fontsize (1,1) double {mustBePositive, mustBeInteger} = 16;
    
    options.rot_theta_zyx (1,3) double = [0 0 0];

    options.disp_mag (1,1) logical = true;
    options.disp_phase (1,1) logical = true;
    options.disp_cb (1,1) logical = true;

    options.disp_title  (1,1) logical = true;
    options.disp_xlabel (1,1) logical = true;
    options.disp_ylabel (1,1) logical = true;
    options.disp_xticks (1,1) logical = true;
    options.disp_yticks (1,1) logical = true;
    
    options.dB_lim (1,2) double  = [-inf inf];
    options.deg_lim (1,2) double = [-180 180];

    options.disp_xaxis_ker_size (1,1) double {mustBePositive, mustBeInteger} = 1;

    options.title_name (1,:) char = [];
    options.title_name_override (1,:) char = [];

    options.fig_pos  (1,2) double = [100, 100];
    options.fig_size (1,2) double {mustBePositive} = [560, 420];
    
    options.disp_theta_phi (:,2) double = [];
    options.disp_theta_phi_markersize (1,1) double {mustBeNonnegative} = 48;
    
end

[P, M] = size(C);
P = floor(sqrt(P))-1;
assert(P - floor(P) == 0, 'Invalid size C');

%Rotate field
if ~isequal(options.rot_theta_zyx, [0, 0, 0])
     C = sh_rot(C, options.rot_theta_zyx, is_real);
end

%Wrap azimuth to be [-pi, pi]
if ~isempty(options.disp_theta_phi)
    options.disp_theta_phi(:, 2) = mod(options.disp_theta_phi(:, 2), 2*pi);
    idx = options.disp_theta_phi(:, 2) > pi;
    options.disp_theta_phi(idx, 2) =  options.disp_theta_phi(idx, 2) - 2 * pi;
end

if strcmp(mode, 'mercator') %Plot mercator layout per function

    phi_list = linspace(-pi, pi, options.resolution(2));
    theta_list = linspace(0, pi, options.resolution(1));
    [phi, theta] = meshgrid(phi_list, theta_list);
    phi = phi(:);
    theta = theta(:);
    Y = sh_val(P, theta, phi, is_real);  %[N x (P+1)^2] 

    h = cell(1, M);
    for m = 1:M

        f = reshape(Y * C(:, m), [options.resolution(1), options.resolution(2)]);

        h{m} = figure;
        h{m}.Position = [options.fig_pos, options.fig_size];
        if options.disp_mag && options.disp_phase
            tiledlayout(2,1);
        elseif options.disp_mag || options.disp_phase
            tiledlayout(1,1);
        else
            error('options.disp_mag or options.disp_phase must be true');
        end

        if options.disp_mag %Plot magnitude response
            nexttile;
            
            h_pc = pcolor(rad2deg(phi_list), 90 - rad2deg(theta_list), mag2db(abs(f)) );
            set(h_pc, 'EdgeColor', 'none');
            set(h_pc, 'FaceColor', options.FaceColor);

            if ~isempty(options.disp_theta_phi)
                hold on;
                scatter(rad2deg(options.disp_theta_phi(:, 2)), 90 - rad2deg(options.disp_theta_phi(:, 1)), options.disp_theta_phi_markersize, 'r', 'linewidth', 2);
            end

            if options.disp_xlabel
                xlabel('Azimuth (Degrees)', 'fontsize', options.fontsize);
            end
            if ~options.disp_xticks
                set(gca,'XTick', []);
            end
            if options.disp_ylabel
                ylabel('Elevation', 'fontsize', options.fontsize);
            end
            if ~options.disp_yticks
                set(gca,'YTick', []);
            end

            if ~isempty(options.hz) %Index at frequency                
                title_str = ['Magnitude Response: ', num2str(options.hz(m)), ' Hz'];
            elseif ~isempty(options.t) %Index at time
                title_str = ['Magnitude Response: ', num2str(options.t(m)), ' Second'];
            else
                if M == 1
                    title_str = ['Magnitude Response'];
                else
                    title_str = ['Magnitude Response: m = ', num2str(m)];
                end                
            end
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

            clim([options.dB_lim(1), options.dB_lim(2)]);
            if options.disp_cb
                h_cb = colorbar;
                ylabel(h_cb, 'dB', 'fontsize', options.fontsize);
            end
        end

        if options.disp_phase %Plot phase response
            nexttile;
            
            h_pc = pcolor(rad2deg(phi_list), 90 - rad2deg(theta_list), rad2deg(angle(f)) );
            set(h_pc, 'EdgeColor', 'none');
            set(h_pc, 'FaceColor', options.FaceColor);

            if ~isempty(options.disp_theta_phi)
                hold on;
                scatter(rad2deg(options.disp_theta_phi(:, 2)), 90 - rad2deg(options.disp_theta_phi(:, 1)), options.disp_theta_phi_markersize, 'w', 'linewidth', 2);
            end

            if options.disp_xlabel
                xlabel('Azimuth (Degrees)', 'fontsize', options.fontsize);
            end
            if ~options.disp_xticks
                set(gca,'XTick', []);
            end
            if options.disp_ylabel
                ylabel('Elevation', 'fontsize', options.fontsize);
            end
            if ~options.disp_yticks
                set(gca,'YTick', []);
            end           
         
            if ~isempty(options.hz) %Index at frequency                
                title_str = ['Magnitude Response: ', num2str(options.hz(m)), ' Hz'];
            elseif ~isempty(options.t) %Index at time
                title_str = ['Magnitude Response: ', num2str(options.t(m)), ' Second'];
            else
                if M == 1
                    title_str = ['Magnitude Response'];
                else
                    title_str = ['Magnitude Response: m = ', num2str(m)];
                end                
            end
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

            clim([options.deg_lim(1), options.deg_lim(2)]);
            if options.disp_cb
                h_cb = colorbar;
                ylabel(h_cb, 'Degrees', 'fontsize', options.fontsize);
            end
            colormap(gca, hsv);

        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode, 'horizontal') %Plot horizontal plane

        phi = linspace(-pi, pi, options.resolution(2))';
        theta = pi/2 * ones(numel(phi), 1);

        Y = sh_val(P, theta, phi, is_real);  %[N x (P+1)^2] 
        f =  Y * C;     %[N x M]
        
        %Optional 1D convolution in time
        if options.disp_xaxis_ker_size > 1
           f = conv2(f, hann(options.disp_xaxis_ker_size)');
           f = f(:, 1:M);
        end

        h = cell(1, 1);
        h{1} = figure;
        h{1}.Position = [options.fig_pos, options.fig_size];
        if options.disp_mag && options.disp_phase
            tiledlayout(2,1);
        elseif options.disp_mag || options.disp_phase
            tiledlayout(1,1);
        else
            error('options.disp_mag or options.disp_phase must be true');
        end

        if options.disp_mag %Plot magnitude response
            nexttile;

            if ~isempty(options.hz) %Index in frequency
                assert(numel(options.hz) == M, 'options.hz size mismatch C');
                h_pc = pcolor(options.hz, rad2deg(phi), mag2db(abs(f)) );
            elseif ~isempty(options.t) %Index in time
                assert(numel(options.t) == M, 'options.t size mismatch C');
                h_pc = pcolor(options.t, rad2deg(phi), mag2db(abs(f)) );
            else
                h_pc = pcolor(1:M, rad2deg(phi), mag2db(abs(f)) );
            end
            set(h_pc, 'EdgeColor', 'none');
            set(h_pc, 'FaceColor', options.FaceColor);

            if options.disp_xlabel
                if ~isempty(options.hz) %Index in frequency
                     xlabel('Frequency (Hz)', 'fontsize', options.fontsize);
                elseif ~isempty(options.t) %Index in time
                     xlabel('Time (Second)', 'fontsize', options.fontsize);
                else                    
                     xlabel('Sample', 'fontsize', options.fontsize);                  
                end
            end
            if ~options.disp_xticks
                set(gca,'XTick', []);
            end
            if options.disp_ylabel
                ylabel('Azimuth (Degrees)', 'fontsize', options.fontsize);
                %ylabel('Azimuth', 'fontsize', options.fontsize);
            end
            if ~options.disp_yticks
                set(gca,'YTick', []);
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

            clim([options.dB_lim(1), options.dB_lim(2)]);
            if options.disp_cb
                h_cb = colorbar;
                ylabel(h_cb, 'dB', 'fontsize', options.fontsize);
            end
            %xscale log;

        end
    
        if options.disp_phase %Plot phase response
            nexttile;

            if ~isempty(options.hz) %Index in frequency
                assert(numel(options.hz) == M, 'options.hz size mismatch C');
                h_pc = pcolor(options.hz, rad2deg(phi), rad2deg(angle(f)) );
            elseif ~isempty(options.t) %Index in time
                assert(numel(options.t) == M, 'options.t size mismatch C');
                h_pc = pcolor(options.t, rad2deg(phi), rad2deg(angle(f)) );
            else
                h_pc = pcolor(1:M, rad2deg(phi), rad2deg(angle(f)) );
            end
            set(h_pc, 'EdgeColor', 'none');
            set(h_pc, 'FaceColor', options.FaceColor);

            if options.disp_xlabel
                if ~isempty(options.hz) %Index in frequency
                     xlabel('Frequency (Hz)', 'fontsize', options.fontsize);
                elseif ~isempty(options.t) %Index in time
                     xlabel('Time (Second)', 'fontsize', options.fontsize);
                else                    
                     xlabel('Sample', 'fontsize', options.fontsize);                  
                end
            end
            if ~options.disp_xticks
                set(gca,'XTick', []);
            end
            if options.disp_ylabel
                ylabel('Azimuth (Degrees)', 'fontsize', options.fontsize);
            end
            if ~options.disp_yticks
                set(gca,'YTick', []);
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

            clim([options.deg_lim(1), options.deg_lim(2)]);
            if options.disp_cb
                h_cb = colorbar;
                ylabel(h_cb, 'Degrees', 'fontsize', options.fontsize);
            end
            %xscale log;

        end

else
    error('Unsupported mode');
end
