function [h, err, h_figs] = ft_fit_sh(C, w, is_real, num_taps, mode, options)
%Regularized least squares fit of FIR filters h to SH coefficients C:
%min_h ||F(w) * (h') - (C.')||^2 + h*Q*h',
%h is real unknowns,
%F is Fourier transform matrix at angular frequency w, 
%C is SH coefficients,
%Q is regularization matrix
%Solution: h = real(F' * F + Q) \ real(F' * (C.'))

%mode = 'identity'  unity Q, penalize squared Euclidean norm of h
%Q = lambda * I

%mode == 'exp'    Complementary Gaussian shaped window for regularization of h 
%Q = lambda * diag(q),     q(n) = 1 - exp( - ( (n-1) - exp_mu )^2 / (2 * exp_std^2) )

%mode == 'circularexp'    Complementary circular Gaussian shaped window for regularization of h 
%Q = lambda * diag(q),     q(n) = c - sum_k exp( - ( (n-1) - exp_mu + k * num_taps )^2 / (2 * exp_std^2) ), k = -inf to inf, c = max(exp_sum)

%mode == 'lognorm'    Complementary log-normal shaped window for regularization of h 
%Q = lambda * diag(q),     q(n) = 1 - exp( - ( log(n-1) - log_mu )^2 / (2 * log_var) )
%where log_mu, log_var are designed from density function's mode and var

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:         [(P + 1)^2 x M] SH coefficients over M frequencies
%w:         [1 x M] Angular frequency (radians / sample) [0 to 2 * pi]
%is_real:   Logical, if true, evaluate real SH
%mode:      String, regularization method {'identity', 'exp', 'circularexp', 'lognorm'}

%options:               Struct

%options.exp_lambda:  Scalar regularization weights for mode = 'exp', 'circularexp'
%options.exp_mu:      Mean delay regularization weights for mode = 'exp', 'circularexp'
%options.exp_std:     Standard deviation delay regularization weights for mode = 'exp', 'circularexp'

%options.lognorm_lambda:  Scalar regularization term weights for mode = 'lognorm'
%options.lognorm_mode:    Mode delay regularization weights for mode = 'lognorm'
%options.lognorm_std:     Standard deviation delay regularization weights for mode ='lognorm'

%options.enable_disp:       Logical, if true, plot filter and fitting error
%options.disp_M_log:        Number of log-frequency sampling points for display
%options.disp_dB_lim:       [1 x 2] dB range [min, max]   
%options.disp_err_dB_lim:   [1 x 2] dB error range [min, max]
%options.Fs:                Sampling rate for display

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output 
%h:         [(P + 1)^2 x num_taps]  Row-matrix filter taps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Resample log-frequency response of spherical piston to
%linear-frequency

% Fs = 48000;
% hz = logspace(log10(20), log10(Fs / 2), 128 + 1);
% w  = hz / Fs * (2 * pi);
% r = 3; % Evaluation distance
% radius = 0.25; % Spherical baffle radius
% C = sh_enc_pist_sphere(30, pi/2, 0, radius, deg2rad(20), hz, r);

% num_taps = 1024;
% is_real = false;
% h_fit_exp = ft_fit_sh(C, w, is_real, num_taps, 'exp', 'exp_mu', (r - radius) / 343 * Fs, 'exp_std', 100, ...
%                   'enable_disp', true, 'disp_dB_lim', [-20, 60], 'disp_err_dB_lim', [-80, 10]);
% h_fit_lognorm = ft_fit_sh(C, w, is_real, num_taps, 'lognorm', 'lognorm_mode', (r - radius) / 343 * Fs, 'lognorm_std', 100, ...
%                   'enable_disp', true, 'disp_dB_lim', [-20, 60], 'disp_err_dB_lim', [-80, 10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Resample log-frequency response of normalized spherical piston to
%linear-frequency

% Fs = 48000;
% hz = logspace(log10(20), log10(Fs / 2), 128 + 1);
% w  = hz / Fs * (2 * pi);
% r = 3; % Evaluation distance
% radius = 0.25; % Spherical baffle radius
% C = sh_enc_pist_sphere(30, pi/2, 0, radius, deg2rad(20), hz, r, 'normalize_piston_axis', true);

% num_taps = 1024;
% is_real = false;
% C_fit = ft_fit_sh(C, w, is_real, num_taps, 'exp', 'exp_std', 150, ...
%                   'enable_disp', true, 'disp_dB_lim', [-60, 6], 'disp_err_dB_lim', [-80, 10]);


arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
    w (1,:) double = 0;
    is_real (1,1) logical = false;
    num_taps (1,1) double {mustBePositive} = 1;
    mode (1,:) char {mustBeMember(mode, {'identity', 'exp', 'circularexp', 'lognorm'})} = 'exp';

    options.exp_lambda (1,1) double {mustBeNonnegative} = 1;
    options.exp_mu (1,1) double = 0;
    options.exp_std (1,1) double {mustBePositive} = 100;

    options.lognorm_lambda (1,1) double {mustBeNonnegative} = 1;
    options.lognorm_mode (1,1)  double  {mustBeNonnegative} = 1;
    options.lognorm_std (1,1) double {mustBePositive} = 50;

    options.enable_disp (1,1) logical = false;
    options.disp_M_log (1,1) {mustBePositive, mustBeInteger} = 2048;
    options.disp_dB_lim (1,2) double = [-inf, inf];
    options.disp_err_dB_lim (1,2) double = [-inf, inf];
    options.Fs (1,1) double {mustBePositive} = 48000;
end

[N_C, M] = size(C);
assert(numel(w) == M, 'C, w size mismatch M');

% Compute Discrete Fourier transform matrix
ndx = 0:(num_taps - 1);
F = exp( -1i * (w(:) * ndx) ); %[M x num_taps]

% Compute regularization matrix Q
if any( strcmp(mode, {'identity', 'exp', 'circularexp', 'lognorm'} ) )

    if strcmp(mode, 'identity')
    
        q = ones(1, num_taps);
        Q = options.exp_lambda * diag(q);
    
    elseif strcmp(mode, 'exp')
    
        q = 1 - exp( - (ndx - options.exp_mu).^2 / (2 * options.exp_std.^2) );
        Q = options.exp_lambda * diag(q);

    elseif strcmp(mode, 'circularexp')
   
        exp_sum = zeros(1, num_taps);
        z_std3 = ceil(3 * options.exp_std / num_taps); % 3 standard deviations
        for i = -z_std3:z_std3
            exp_sum = exp_sum + exp( - (ndx - options.exp_mu + i * num_taps).^2 / (2 * options.exp_std.^2) );
        end
        q = max(exp_sum) - exp_sum;
        Q = options.exp_lambda * diag(q);

    elseif strcmp(mode, 'lognorm')
    
        %(u - 1) * u^3 - var_x^2 / mode_x^2 = 0, where u = exp(var)
        u = roots([1, -1, 0, 0, -options.lognorm_std^2 / options.lognorm_mode^2]);
        u_real_roots = u(abs(imag(u)) <= 1e-8 );
        u_pos_real = u_real_roots(u_real_roots > 0);
        u_pos_real = u_pos_real(1);

        if ~isempty(u_pos_real)
            log_var = log(u_pos_real);
            log_mu  = log(options.lognorm_mode) + log_var;

            q = 1 - exp( - (log(ndx) - log_mu).^2 / (2 * log_var) );
            Q = options.lognorm_lambda * diag(q);
        else
            error('Invalid lognorm_mode or lognorm_std');
        end

    else
        error('Unsupported mode');
    end

    % Compute real-valued least squares solution 
    h = (real(F' * F + Q) \ real(F' * C.')).'; % LHS is Hermitian

else
    error('Unsupported mode');
end

% Compute error
err = norm(F * (h.') - (C.'));

% Plotting
h_figs = [];
if options.enable_disp && coder.target("MATLAB")

    hz_w = w / (2 * pi) * options.Fs;
    hz_log = logspace(log10(20), log10(options.Fs / 2), options.disp_M_log);

    H_hz_log = zeros([N_C, options.disp_M_log]);
    for n = 1:N_C
        H_hz_log(n, :) = freqz(h(n, :), 1, hz_log, options.Fs);
    end

    H_hz_w = zeros([N_C, M]);
    for n = 1:N_C
        H_hz_w(n, :) = freqz(h(n, :), 1, hz_w, options.Fs);
    end

    % Display filters
    fontsize = 16;
    h_filters = figure;
    yyaxis left;
    plot(ndx(:), h.', '-', 'Color', [0, 0, 1, max(0.01, 10/N_C)], 'linewidth', 1.5);
    ylabel('Filter Taps', 'fontsize', fontsize);
    yyaxis right;
    plot(ndx(:), diag(Q), 'linewidth', 1.5);
    ylabel('Regu. Var.', 'fontsize', fontsize);
    xlabel('Time samples', 'fontsize', fontsize);
    title('Fitted Filters', 'fontsize', fontsize + 1);
    grid on; axis tight;
    set(gca, 'fontsize', fontsize - 1);

    % Display target and fit
    h_fig_tgt = sh_plt(C, 'horizontal', is_real, 'title_name', 'Target', 'hz', hz_w, 'dB_lim', options.disp_dB_lim);
    h_fig_fit = sh_plt(H_hz_log, 'horizontal', is_real, 'title_name', 'Fit', 'hz', hz_log, 'dB_lim', options.disp_dB_lim);
    h_fig_res = sh_plt(H_hz_w - C, 'horizontal', is_real, 'title_name', 'Residual: Fit - Target', 'hz', hz_w, 'dB_lim', options.disp_err_dB_lim);

    h_figs = { h_filters, h_fig_tgt{1}, h_fig_fit{1}, h_fig_res{1} };
end