function [D, h_figs] = sh_translate(C, omega, is_real, v, options)
% Translate SH expansion C's origin to coordinate v. Refit to higher max-order expansion D.

% Restrict translated field to delay-only modifications under far-field assumptions. 
% SH C encodes planewave frequency responses from incident direction u.
% SH D encodes planewave frequency responses delayed w.r.t. the new origin under projection of 
% translation vector v onto normal direction u.

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:        [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%omega:    [1 x M] Angular frequency (2 * pi * hz)
%is_real:  Logical, if true, evaluate real SH
%v:        [1 x 3] Cartesian coordinate of new origin (meters)

%options.c:     Wavelength at 1 Hz in meters / second
%options.R:     Radius of scatterer or geometry

%options.oversample: Oversample factor (1 = none)
%options.max_odr_fit:   Maximum order SH fit

%options.enable_disp:   Logical, if true, display fields
%options.disp_dB_lim:   [1 x 2] Display dB range [min, max]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%D:        [(P_E + 1)^2 x M] SH coefficients where P_E = min( options.max_odr_fit, max(P, floor(max(kappa) * (options.R + norm(v))) ) )
%                            for wavenumber kappa = omega / options.c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate SH spherical piston normalized on-axis and translate the expansion

% Fs = 8000;
% is_real = false;
% N_taps = 128;
% w = linspace(0, 2 * pi, N_taps + 1); w = w(1:end-1);
% t = (0:(N_taps-1))/Fs;
% omega = w * Fs;
% freq = omega / (2 * pi);

% C = sh_enc_pist_sphere(30, pi/2, 0, 0.25, deg2rad(20), freq, 4, 'normalize_piston_axis', true);
% D = sh_translate(C, omega, is_real, [- 0.5/1000 * 343, 0, 0], 'enable_disp', true, 'disp_dB_lim', [-20, 6]); % 1/2 ms in -x direction
% 
% sh_plt(C, 'cardinal', is_real, 'hz', freq, 'title_name', 'Original');
% sh_plt(D, 'cardinal', is_real, 'hz', freq, 'title_name', 'Translated');

% C_td = ifft(C, [], 2, 'symmetric');
% sh_plt(C_td, 'cardinal', is_real, 't', t, 'title_name', 'Original');
% D_td = ifft(D, [], 2, 'symmetric');
% sh_plt(D_td, 'cardinal', is_real, 't', t, 'title_name', 'Translated');

arguments
    C (:,:) double = complex(0);
    omega (1,:) double = 0;
    is_real (1,1) logical = false;
    v (1,3) double = [0 0 0];

    options.c (1,1) double {mustBePositive} = 343;  % Sound in air under normative temp and humidity
    options.R (1,1) double {mustBePositive} = 0.01; % Radius of bounding sphere
    options.oversample (1,1) double {mustBeInteger, mustBePositive} = 4;
    options.max_odr_fit (1,1) double {mustBeInteger, mustBeNonnegative} = 60;

    options.enable_disp (1,1) logical = false;
    options.disp_dB_lim (1,2) double = [-inf, inf];
end

[N_C, M] = size(C);
max_odr = sqrt(N_C) - 1;
assert(max_odr - floor(max_odr) == 0, 'Invalid size C');

kappa = omega / options.c; % Wavenumber
P_max = floor(max(kappa) * (options.R + norm(v))); % Theoretical maximum expansion order
P_E = min( options.max_odr_fit, max(P_max, max_odr) );

[theta, phi, U] = sh_grd_fib( options.oversample * (P_E + 1)^2 ); %Over-sample

d_meters = U * v'; % [N_C x 1] Signed distance of v projected onto unit directions U in meters
d_sec = d_meters / options.c; % Signed delay in seconds
    
f = sh_dec(C, theta, phi, is_real); % [N_C x M]
f_del = f .* exp(1i * d_sec * omega); % Delayed frequency responses

D = sh_fit_svd(f_del, theta, phi, P_E, is_real);

% Plotting
h_figs = [];
if options.enable_disp && coder.target("MATLAB")

    freq = omega / (2 * pi);
    
    h_orig = sh_plt(C(:, freq > 0), 'horizontal', is_real, 'title_name', 'Original', 'dB_lim', options.disp_dB_lim, 'hz', freq(freq > 0));
    h_trans = sh_plt(D(:, freq > 0), 'horizontal', is_real, 'title_name', 'Translated', 'dB_lim', options.disp_dB_lim, 'hz', freq(freq > 0));
    
    h_figs = {h_orig{1}, h_trans{1}};
end