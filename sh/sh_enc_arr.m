function [C, center, h_figs] = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, options)
% Encode planewave to analytic transducer (e.g. microphone, loudspeaker, antenna) array 
% frequency responses in spherical harmonic domain. SH expansion centers across transducers
% are at the common origin coordinates [0, 0, 0].

% If options.noncentered is true, SH expansion centers are moved to each transducuer's coordinates and compensated later 
% via output center coordinates when decoding. Uniform directivity arrays are therefore encoded by just their 0th order expansions
% as direction-dependent delays due to coordinate translations become constant.

% Author:    Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%max_odr:       Max SH order
%mode:          String, array element type and topology {'uni_line', 'uni_circ', 'uni_sphere', 'uni_custom', 'spis_circ', 'spis_sphere', 'spis_custom'}
%                   'uni_line':     Uniform (omni) directivity, line array on x=z=0 line
%                                   equal spaced with first element on -y axis
%                   'uni_circ':     Uniform (omni) directivity, circular topology on z=0 plane,
%                                   equal angle element spacing with first element on +x axis
%                   'uni_sphere':   Uniform (omni) directivity, Fibonacci sphere points
%                   'uni_custom':   Uniform (omni) directivity, specify cartesian coordinates in options.custom_cartesian, ignore radius
%                   'spis_circ':    Piston on spherical baffle directivity, circular topology on z=0 plane
%                   'spis_sphere':  Piston on spherical baffle directivity, Fibonacci sphere points
%                   'spis_custom':  Piston on spherical baffle directivity, specify spherical coordinates options.custom_theta_phi

%                   'fit_rot_custom':   Rotate (intrinsic Tait-Bryan) fitted expansions, specify expansion in options.custom_C, 
%                                       spherical coordinates options.custom_theta_phi, and rotation angles options.custom_rot_zyx.
%                                       Expansion center cshould be at rotationally invariant coordinate (e.g. center of a spherical baffle).

%N:             Number of array elements
%radius:        Array radius in meters
%delay_dist:    Delay distance to apply to all frequency responses in meters
%freq:          [1 x M] Frequency in Hz

%options.c:             Velocity of propagation meters / second

%options.custom_cartesian:  [3 x N] custom cartesian coordinates for mode = 'uni_custom'
%options.custom_theta_phi:  [2 x N] custom spherical coordinates (radians) for mode = 'spis_custom'
%options.custom_rot_zyx:    [3 x N] custom zyx intrinsic Tait-Bryan rotation angles (radians) for mode = 'fit_rot_custom'

%options.noncentered:   Logical, if true, use non-centered expansion center coordinates

%options.spis_arc_len:  Spherical piston's arc length in meters for mode = 'spis_circ'

%options.enable_disp:   Logical, if true, display fields
%options.disp_dB_lim:   [1 x 2] Display dB range [min, max]
%options.disp_mag:      Logical, if true, display magnitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:             [(max_odr + 1)^2 x N x M] SH expansions
%center:        [] or [N x 3] transducer dependent expansion coordinates if options.noncentered is true
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
% C = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, 'enable_disp', true, 'disp_dB_lim', [-3, 3]);
% 
% t = 0:(1/Fs):((M-1)/Fs);
% for n = 1:N
%     sh_plt(C(:, n, :), 'cardinal', false, 't', t, 'hz', freq, 'title_name', ['Trans. ', num2str(n)]);
% end

% %Time-domain from DC to Nyquist only
% for n = 1:N
%     sh_plt(C(:, n, 1:(floor(M/2)+1)), 'cardinal', false, 't', t, 'hz', freq, 'title_name', ['Trans. ', num2str(n)]);
% end

% Compare to non-center expansion with max_odr 0
% [C_noncenter, center] = sh_enc_arr(0, mode, N, radius, delay_dist, freq,  'noncentered', true, 'enable_disp', true, 'disp_dB_lim', [-3, 3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate 4 element piston directivity on spherical baffle array responses

% Fs = 24000;
% M = 256;
% N = 4;
% radius = 0.01;
% delay_dist = 1;
% max_odr = 36;
% mode = 'spis_circ';
% freq = linspace(0, Fs, M + 1); freq = freq(1:end-1);
% spis_arc_len = 0.0035; % 3.5 mm
% 
% C = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, 'spis_arc_len', spis_arc_len, 'enable_disp', true, 'disp_dB_lim', [-24, 6]);
% 
% t = 0:(1/Fs):((M-1)/Fs);
% for n = 1:N
%     sh_plt(C(:, n, :), 'cardinal', false, 't', t, 'hz', freq, 'title_name', ['Trans. ', num2str(n)]);
% end

% %Compare to non-center expansions for larger spherical baffle and piston size
% radius = 0.1;
% max_odr = 40;
% spis_arc_len = 0.035; % 3.5 cm
% freq = logspace(log10(20), log10(Fs), M + 1); freq = freq(1:end-1);
% C = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, 'spis_arc_len', spis_arc_len, 'enable_disp', true, 'disp_dB_lim', [-24, 12]);
% [C_noncenter, center] = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, 'noncentered', true, 'spis_arc_len', spis_arc_len, 'enable_disp', true, 'disp_dB_lim', [-24, 12]);

arguments
    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;
    mode (1,:) char {mustBeMember(mode, {'uni_line', 'uni_circ', 'uni_sphere', 'uni_custom', 'spis_circ', 'spis_sphere', 'spis_custom', 'fit_rot_custom'} )} = 'uni_circ'
    N (1,1) double {mustBePositive, mustBeInteger} = 3;
    radius (1,1) double {mustBeNonnegative} = 0.1;
    delay_dist (1,1) double {mustBeNonnegative} = 1;
    freq (1,:) double {mustBeNonnegative} = 0;
    
    options.c (1,1) double {mustBePositive} = 343;

    options.custom_cartesian (3,:) double = [];
    options.custom_theta_phi (2,:) double = [];
    options.custom_rot_zyx (3,:) double = [];
    options.custom_C (:,:) double {coder.mustBeComplex} = [0];
 
    options.noncentered (1,1) logical = false;

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

        assert(size(options.custom_cartesian, 2) == N, 'options.custom_cartesian size mismatch N');
        V = options.custom_cartesian;

    else
        error('Unknown mode');
    end

    d_sec = delay_dist / options.c;  
    C_origin = repmat( sh_enc_uni(0), [1, M]); % uniform directivity at origin
    C_origin = bsxfun(@times, C_origin, exp(-1i * d_sec * omega) );   % Add bulk-delay to make responses causal    

    if options.noncentered % Omit SH translation

        for n = 1:N
            C(:, n, :) = sh_resize(C_origin, max_odr);
        end    
        center = V';

    else

        % Translate expansion origin to transducer coordinates
        for n = 1:N
            C(:, n, :) = sh_resize(sh_translate(C_origin, omega, false, V(:, n)', 'max_odr_fit', max_odr, 'R', radius), max_odr);    
        end
        center = [];
        
    end


elseif any( strcmp(mode, {'spis_circ', 'spis_sphere', 'spis_custom'}) ) % Piston on spherical baffle directivity
   
    theta_breadth = options.spis_arc_len / radius;
    assert(theta_breadth <= 2 * pi, 'arc_len exceeds circumference of spherical baffle');

    d_sec = delay_dist / options.c; 
    C_top = sh_enc_pist_sphere(max_odr, 0, 0, radius, theta_breadth, freq, delay_dist + radius, 'normalize_piston_axis', true);  

    if strcmp(mode, 'spis_circ') 

        phi_list = deg2rad(0:(360/N):(360 - (360/N)) );
        theta_list = pi/2 * ones(1, N);
        radius_list = radius * ones(1,N);

        % Equally distribute coordinates on z=0 plane circle
        [V(1,:), V(2,:), V(3,:)] = sph2cart( phi_list, pi/2 - theta_list, radius_list );        

    elseif strcmp(mode, 'spis_sphere')
        
        % Fibonacci spherical grid
        [theta_list, phi_list, V_transposed] = sh_grd_fib(N);
        theta_list = theta_list';
        phi_list = phi_list';
        V = radius * V_transposed';

    elseif strcmp(mode, 'spis_custom')

        assert(size(options.custom_theta_phi, 2) == N, 'options.custom_theta_phi size mismatch N');
        
        theta_list = options.custom_theta_phi(1,:);
        phi_list = options.custom_theta_phi(2,:);
        
        [V(1,:), V(2,:), V(3,:)] = sph2cart( phi_list, pi/2 - theta_list, radius * ones(1,N) );  

    else
        error('Unknown mode')
    end

    if options.noncentered % sh_enc_pist_sphere expansion center is already at origin, translate in inverse direction

        C_top = sh_resize(sh_translate(C_top, omega, false, [0, 0, -radius], 'max_odr_fit', max_odr, 'R', radius), max_odr);
        center = V';

    else

        center = [];

    end
    C_top = bsxfun(@times, C_top, exp(-1i * d_sec * omega) );   % Add bulk-delay to make responses causal 

    % Rotate expansion onto transducer directions
    for n = 1:N
        C(:, n, :) = sh_rot(C_top, [phi_list(n), theta_list(n), 0], false);
    end 

elseif any( strcmp(mode, {'fit_rot_custom'}) ) % Fitted expansions

    if strcmp(mode, 'fit_rot_custom') % Rotate fitted expansions
        
        assert(size(options.custom_C, 2) == M, 'options.custom_C size mismatch M');
        assert(size(options.custom_theta_phi, 2) == N, 'options.custom_theta_phi size mismatch N');
        assert(size(options.custom_rot_zyx, 2) == N, 'options.custom_rot_zyx size mismatch N');

        theta_list = options.custom_theta_phi(1,:);
        phi_list = options.custom_theta_phi(2,:);        
        [V(1,:), V(2,:), V(3,:)] = sph2cart( phi_list, pi/2 - theta_list, radius * ones(1,N) );  

        C_res = sh_resize(options.custom_C, max_odr);

        if options.noncentered % expansion center is already at origin, translate in inverse direction

            C_res = sh_resize(sh_translate(C_res, omega, false, [0, 0, -radius], 'max_odr_fit', max_odr, 'R', radius), max_odr);
            center = V';

        else

            center = [];

        end

        d_sec = delay_dist / options.c; 
        C_res = bsxfun(@times, C_res, exp(-1i * d_sec * omega) );   % Add bulk-delay to make responses causal 

        % Duplicate and rotate expansions
        C = repmat(reshape(C_res, [size(C_res, 1), 1, M] ), [1, N, 1]);        
        for n = 1:N
            C(:, n, :) = sh_rot(options.custom_C, options.custom_rot_zyx(:, n)', false);
        end

    else
        error('Unknown mode');
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

    if contains(mode, 'spis_')
        [sx, sy, sz] = sphere(200);
        surf(sx * radius, sy * radius, sz * radius, 'FaceAlpha', 0.2, 'FaceColor', 'interp', 'EdgeColor', 'none'); 
        for n = 1:N
            V_unit_n = V(:,n)' / radius;
            radius_n = sin(theta_breadth) * radius;
            plot_great_circle(V_unit_n, radius, radius_n); hold on;
        end
    end

    xlabel('X-axis (meters)', 'fontsize', fontsize);
    ylabel('Y-axis (meters)', 'fontsize', fontsize);
    zlabel('Z-axis (meters)', 'fontsize', fontsize);
    title('Transducer 3D Coordinates', 'fontsize', fontsize + 1);
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
        str_title_n = ['Trans. ', num2str(n), ' [', num2str(V(1,n),3), ' ,', num2str(V(2,n),3), ', ', num2str(V(3,n),3), ']'];
        
        idx_pos_freq = find(freq > 0);
        if isempty(center)
            center_wavnum = [];
        else
            center_wavnum = [repmat(center(n, :), [M, 1]), omega(:) / options.c];            
            center_wavnum = center_wavnum(idx_pos_freq, :);
        end
        h_n = sh_plt(reshape(C(:, n, idx_pos_freq), [N_C, numel(idx_pos_freq)]), 'horizontal', false, 'center_wavnum', center_wavnum, ...
            'title_name', str_title_n, 'dB_lim', options.disp_dB_lim, 'hz', freq(idx_pos_freq), 'disp_mag', options.disp_mag);
        h_figs{end+1} = h_n{1};
    end    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_great_circle(v, r_sphere, r_circle)
% Plot great circle in direction of unit-vector v, of radius r_ricle in sphere of radius r_sphere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% v:            [1 x 3] Unit-vector direction of great circle
% r_sphere:     Radius of sphere
% r_circle:     Radius of great circle

arguments
    v (1,3) double = [1 0 0];
    r_sphere (1,1) double {mustBeNonnegative} = 1;
    r_circle (1,1) double {mustBeNonnegative} = 0.1;
end

assert(r_circle <= r_sphere, 'r_circle must be <= r_sphere');

theta = linspace(0, 2 * pi, 100);
U = null(v);  % Unit-vectors orthogonal to v
d = sqrt(r_sphere^2 - r_circle^2);
coords = bsxfun(@plus, r_circle * (U(:,1) * cos(theta) + U(:,2) * sin(theta)), d * v(:));

plot3( coords(1,:), coords(2,:), coords(3,:), 'r-', 'linewidth', 1.25 );

