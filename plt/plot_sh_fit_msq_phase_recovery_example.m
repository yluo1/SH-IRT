function plot_sh_fit_msq_phase_recovery_example(max_odr, options)
% Recover phase response from magnitude squared response of random SH field

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%max_odr:           Max-order expansion in reference field   

%options:               struct
%options.rseed:         Randomized seed
%options.oversample:    Oversample spherical coordinates factor (1 = none)
%options.dB_lim:        [1 x 2] Display dB range [min, max]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Sample phase recovery from of random field from its magnitude squared

% plot_sh_fit_msq_phase_recovery_example(3);
% plot_sh_fit_msq_phase_recovery_example(4);
% plot_sh_fit_msq_phase_recovery_example(5);

arguments
    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 3;
    options.rseed (1,1) double {mustBeNonnegative, mustBeInteger} = 5342; 
    options.oversample (1,1) double {mustBePositive, mustBeInteger} = 3;
    options.dB_lim (1,2) double = [-40, 20];
end

rng(options.rseed);
is_real = false;

C_ref = sh_rand(max_odr, 1, is_real); % Random field
%C_ref = sh_enc_pist_sphere(max_odr, pi/2, 0, 0.25, deg2rad(20), 1000, 4); % Spherical piston

[theta, phi] = sh_grd_fib(options.oversample * (2 * max_odr + 1)^2); 
X = sh_dec(C_ref, theta, phi, is_real);
X_MS = abs(X).^2; % Magnitude Squared target

% Magnitude squared fit
[D_fit, C_fit] = sh_fit_msq(X_MS, theta, phi, 2 * max_odr, is_real, 'MS', 'C0', sh_rand(max_odr, options.oversample * (2 * max_odr + 1)^2, is_real));

% Normalize phase at theta = pi/2, phi = 0 to 0 radians
C_ref = C_ref * exp(-1i * angle(sh_dec(C_ref, pi/2, 0, is_real)));
C_fit = C_fit * exp(-1i * angle(sh_dec(C_fit, pi/2, 0, is_real)));
C_fit_conj = sh_conj(C_fit); % Allow conjugation as alternative candidate

% Compute error
err_fit = err_NSHMSQ(C_ref, C_fit)
err_fit_conj = err_NSHMSQ(C_ref, C_fit_conj)

if err_fit < err_fit_conj
    err = err_fit;
    C_fit_best = C_fit;
else
    err = err_fit_conj;
    C_fit_best = C_fit_conj;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err
h_ref = sh_plt(C_ref, 'mercator', is_real, 'title_name', 'Original', 'dB_lim', options.dB_lim);
h_fit = sh_plt(C_fit_best, 'mercator', is_real, 'title_name', ['Fit NMSQE ', num2str(err)], 'dB_lim', options.dB_lim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_dir = 'figs/figs_sh/fit';
if ~isfolder(out_dir)
    mkdir(out_dir);
end

exportgraphics(h_ref{1}, fullfile(out_dir,  ['msq_ms_phaserecover_ref_', num2str(max_odr), '.png'] ));
exportgraphics(h_fit{1}, fullfile(out_dir, ['msq_ms_phaserecover_fit_', num2str(max_odr), '.png'] ));
