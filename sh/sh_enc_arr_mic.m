function [C, h_figs] = sh_enc_arr_mic(max_odr, mode, N, radius, delay_dist, freq, options)
% Encode planewave to analytic microphone array frequency responses in spherical harmonic domain

% Author:    Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%max_odr:       Max SH order
%mode:          String, array element type and topology {'uni_line', 'uni_circ', 'uni_sphere', 'uni_custom', 'pist_sphere'}
%                   'uni_line':     Uniform (omni) directivity, line array on x=z=0 line
%                                   equal spaced with first element on -y axis
%                   'uni_circ':     Uniform (omni) directivity, circular topology on z=0 plane,
%                                   equal angle element spacing with first element on +x axis
%                   'uni_sphere':   Uniform (omni) directivity, Fibonacci sphere points
%                   'uni_custom':   Uniform (omni) directivity, use options.custom_coords, ignore radius
%                   'spis_circ':    Piston on spherical baffle directivity, circular topology on z=0 plane

%N:             Number of array elements
%radius:        Array radius in meters
%delay_dist:    Delay distance to apply to all frequency responses in meters
%freq:          [1 x M] Frequency in Hz

%options.c:             Velocity of sound meters / second
%options.V_custom:      [3 x N] custom coordinates for mode = 'uni_custom'

%options.spis_arc_len:  Spherical piston's arc length in meters for mode = 'spis_circ'

%options.enable_disp:   Logical, if true, display fields
%options.disp_dB_lim:   [1 x 2] Display dB range [min, max]
%options.disp_mag:      Logical, if true, display magnitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:             [(max_odr + 1)^2 x N x M] SH expansions
%h_figs:        [1 x (N + 1)] cell array of figure handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate 3 element uniform directivity array responses

% Fs = 16000;
% M = 256;
% N = 3;
% radius = 0.1;
% delay_dist = 1;
% max_odr = 30;
% %mode = 'uni_line';
% mode = 'uni_circ';
% %mode = 'uni_sphere';
% freq = linspace(0, Fs, M + 1); freq = freq(1:end-1);
% 
% C = sh_enc_arr_mic(max_odr, mode, N, radius, delay_dist, freq, 'enable_disp', true, 'disp_dB_lim', [-3, 3]);
% 
% t = 0:(1/Fs):((M-1)/Fs);
% for n = 1:N
%     sh_plt(C(:, n, :), 'cardinal', false, 't', t, 'hz', freq, 'title_name', ['Mic ', num2str(n)]);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate 4 element piston directivity on spherical baffle array responses

% Fs = 16000;
% M = 256;
% N = 4;
% radius = 0.01;
% delay_dist = 1;
% max_odr = 36;
% mode = 'spis_circ';
% freq = linspace(0, Fs, M + 1); freq = freq(1:end-1);
% spis_arc_len = 0.0035; % 3.5 mm
% 
% C = sh_enc_arr_mic(max_odr, mode, N, radius, delay_dist, freq, 'spis_arc_len', spis_arc_len, 'enable_disp', true, 'disp_dB_lim', [-6, 6]);
% 
% t = 0:(1/Fs):((M-1)/Fs);
% for n = 1:N
%     sh_plt(C(:, n, :), 'cardinal', false, 't', t, 'hz', freq, 'title_name', ['Mic ', num2str(n)]);
% end

arguments
    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;
    mode (1,:) char {mustBeMember(mode, {'uni_line', 'uni_circ', 'uni_sphere', 'uni_custom', 'spis_circ'} )} = 'uni_circ'
    N (1,1) double {mustBePositive, mustBeInteger} = 3;
    radius (1,1) double {mustBeNonnegative} = 0.1;
    delay_dist (1,1) double {mustBeNonnegative} = 1;
    freq (1,:) double {mustBeNonnegative} = 0;
    
    options.c (1,1) double {mustBePositive} = 343;
    options.custom_coords (3,:) double = [];

    options.spis_arc_len (1,1) double {mustBePositive} = 0.0035;

    options.enable_disp (1,1) logical = false;
    options.disp_dB_lim (1,2) double = [-inf, inf];
    options.disp_mag (1,1) logical = true;
end

M = numel(freq);
N_C = (max_odr + 1)^2;

V = zeros(3, N); % Array element coordinates
C = complex(zeros([N_C, N, M]));

omega = 2 * pi * freq;

if any(strcmp(mode, {'uni_line', 'uni_circ', 'uni_sphere', 'uni_custom'})) % Uniform directivity

    if strcmp(mode, 'uni_line')
        
        % Equi-space coordinates along x=z=0 line
        V(2,:) = linspace(-radius, radius, N);

    elseif strcmp(mode, 'uni_circ')

        % Equally distribute coordinates on z=0 plane circle
        [V(1,:), V(2,:), V(3,:)] = sph2cart( deg2rad(0:(360/N):(360 - (360/N)) ), zeros(1,N), radius * ones(1,N) );
    
    elseif strcmp(mode, 'uni_sphere')
        
        % Fibonacci spherical grid
        [~, ~, V_transposed] = sh_grd_fib(N);
        V = radius * V_transposed';

    elseif strcmp(mode, 'uni_custom')

        assert(size(options.custom_coords, 2) == N, 'options.custom_coords size mismatch N');
        V = options.custom_coords;

    else
        error('Unknown mode');
    end

    d_sec = delay_dist / options.c;  
    C_origin = repmat( sh_enc_uni(0), [1, M]); % uniform directivity at origin

    % Translate expansion origin to microphone coordinates
    for n = 1:N
        C(:, n, :) = sh_resize(sh_translate(C_origin, omega, false, V(:, n)', 'max_odr_fit', max_odr, 'R', radius), max_odr);    
        C(:, n, :) = bsxfun(@times, squeeze(C(:, n, :)), exp(-1i * d_sec * omega) );   % Add bulk-delay to make responses causal    
    end

elseif any( strcmp(mode, {'spis_circ'}) ) % Piston on spherical baffle directivity
   
    theta_breadth = options.spis_arc_len / radius;
    assert(theta_breadth <= 2 * pi, 'arc_len exceeds circumference of spherical baffle');

    d_sec = delay_dist / options.c; 
    C_top = sh_enc_pist_sphere(max_odr, 0, 0, radius, theta_breadth, freq, delay_dist + radius, 'normalize_piston_axis', true);    
    C_top = bsxfun(@times, C_top, exp(-1i * d_sec * omega) );   % Add bulk-delay to make responses causal 

    if strcmp(mode, 'spis_circ') 

        phi_list = deg2rad(0:(360/N):(360 - (360/N)) );
        theta_list = pi/2 * ones(1, N);
        radius_list = radius * ones(1,N);

        % Equally distribute coordinates on z=0 plane circle
        [V(1,:), V(2,:), V(3,:)] = sph2cart( phi_list, pi/2 - theta_list, radius_list );        

    else
        error('Unknown mode')
    end

    % Rotate expansion
    for n = 1:N
        C(:, n, :) = sh_rot(C_top, [phi_list(n), theta_list(n), 0], false);
    end

else
    error('Unknown mode');
end

% Plotting
h_figs = [];
if options.enable_disp && coder.target("MATLAB")
    h_figs = {};
    
    % Plot coordinate map
    h_map = figure;
    h_figs{1} = h_map;

    fontsize = 16;

    scatter3(V(1,:), V(2,:), V(3,:), 48, 'k', 'filled', 'linewidth', 2); hold on;
    constantplane("z", 0, 'FaceColor', 'k', 'FaceAlpha', 0.2);
    constantplane("y", 0, 'FaceColor', 'r', 'FaceAlpha', 0.2);
    constantplane("x", 0, 'FaceColor', 'b', 'FaceAlpha', 0.2);

    xlabel('X-axis', 'fontsize', fontsize);
    ylabel('Y-axis', 'fontsize', fontsize);
    zlabel('Z-axis', 'fontsize', fontsize);
    title('Microphone 3D Coordinates', 'fontsize', fontsize + 1);
    xlim(1.25 * [-radius, radius]);
    ylim(1.25 * [-radius, radius]);
    zlim(1.25 * [-radius, radius]);
    set(gca, 'fontsize', fontsize - 1);
    for n = 1:N
        text(V(1,n), V(2,n), V(3,n), num2str(n), 'VerticalAlignment', 'bottom', 'fontsize', fontsize - 1);
    end
    axis square;
    view(60, 30);

    % Plot horizontal projections
    for n = 1:N
        str_title_n = ['Mic. ', num2str(n), ' [', num2str(V(1,n),3), ' ,', num2str(V(2,n),3), ', ', num2str(V(3,n),3), ']'];        
        h_n = sh_plt(C(:, n, freq > 0), 'horizontal', false, 'title_name', str_title_n, 'dB_lim', options.disp_dB_lim, 'hz', freq(freq > 0), 'disp_mag', options.disp_mag);
        h_figs{end+1} = h_n{1};
    end    

end