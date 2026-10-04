function [C] = sh_filter_val(h, w)
%Compute frequency response of FIR fitted spherical harmonic expansions over at 
%angular frequency w

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%h:         [(P + 1)^2 x num_taps]  Row-matrix filter taps
%w:         [1 x M] Angular frequency (radians / sample) [0 to 2 * pi]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:         [(P + 1)^2 x M] SH coefficients over M frequencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Resample log-frequency response of spherical piston to
%linear-frequency

% Fs = 48000;
% hz = logspace(log10(10), log10(Fs / 2), 256 + 1);
% w  = hz / Fs * (2 * pi);
% r = 3; % Evaluation distance
% radius = 0.25; % Spherical baffle radius
% C = sh_enc_pist_sphere(30, pi/2, 0, radius, deg2rad(20), hz, r);

% is_real = false;
% num_taps = 1024;
% h_fit_exp = ft_fit_sh(C, w, is_real, num_taps, 'exp', 'exp_mu', (r - radius) / 343 * Fs, 'exp_std', 100, ...
%                   'enable_disp', false, 'disp_dB_lim', [-20, 60], 'disp_err_dB_lim', [-80, 10]);

% w_lin = linspace(0, pi, 1024);
% C_lin = sh_filter_val(h_fit_exp, w_lin);
% sh_plt(C, 'horizontal', is_real,  'hz', w / (2 * pi) * Fs, 'dB_lim', [-20, 60], 'title_name', 'Original');
% sh_plt(C_lin, 'horizontal', is_real,  'hz', w_lin / (2 * pi) * Fs, 'dB_lim', [-20, 60], 'title_name', 'Resampled');

arguments
    h (:,:) double {coder.mustBeComplex} = complex(0);
    w (1,:) double = 0;
end

[N_C, num_taps] = size(h);
M = numel(w);
C = zeros([N_C, M]);
for n = 1:N_C
    C_aug = freqz(h(n, :), 1, [0, w]);
    C(n, :) = C_aug(2:end);
end