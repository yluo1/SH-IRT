function [C, t, h_fig] = sh_rand_pp(max_odr, aed, is_real, options)
%Generate Poisson process of random spherical harmonic projections over the spherical coordinates

%Poisson process's inter-arrival times distribute over absolute echo density profile
%in the following reference

%Huang, P. and Abel, J.S., 2007, October. 
%Aspects of reverberation echo density. 
%Audio Engineering Society Convention 123. Audio Engineering Society.

%Randomized projections are uniformly sampled over spherical coordinates unless 
%sampled from non-empty options.C_pdf (SH expansion of probability density function)

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%max_odr:       Max SH order 
%aed:           [1 x M] absolute echo density (# echos / second) profile from 1:M samples
%is_real:       Logical, if true, C is real, otherwise, C is complex, non-negative

%options:       struct
%options.Fs:                    Sampling rate
%options.sinc_window_halflen:   Half-length of Lanczos kernel for sinc interpolation of pulses
%options.pulse_sample_width:    Pulse with in samples

%options.randomize_phase:       Logical, if true, randomize amplitude sign of pulse
%options.normalize_energy:      Logical, if true, normalize to constant energy over time 

%options.direct_unity_first_pulse:    Logical, if true, first pulse is unity at theta = pi/2, phi = 0

%options.C_pdf:                 [(max_odr_pdf + 1)^2 x 1] (independent over M samples)
%                               [(max_odr_pdf + 1)^2 x M] (dependent over M samples)
%                               [] for uniform distribution
%                               SH expansion of spherical probability density function of a projection direction

%options.enable_disp:           Logical, if true, plot figure
%options.disp_dB_lim:           [1 x 2] Display dB [min, max]
%options.disp_xaxis_ker_size:   Kernel size for plotting reflections, positive integer
%options.disp_theta_phi:        [1 x 2] Decode at [co-latitude, azimuth] radians for display

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:     [(max_odr + 1)^2 x M] SH coefficients
%t:     [1 x M] time (sec)
%h_fig: Cell array of figure handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate sample absolute echo density curve and realize a sample path over time
% 
% rng(239);
% Fs = 48000;
% duration = 1;
% M = Fs * duration; %Number of samples
% max_odr = 10;
% aed = logspace(log10(10/Fs), log10(2), M) * Fs;
% is_real = true;

% [C, t] = sh_rand_pp(max_odr, aed, is_real, 'Fs', Fs);
% y0 = real(sh_dec(C(1,:), pi/2, 0)); %Decode 0th order SH
% 
% fontsize = 16;
% figure; yyaxis left;  plot(t, y0, 'linewidth', 1.25);  
% grid on; axis tight; ylabel('Amplitude', 'fontsize', fontsize);
% yyaxis right;  plot(t, aed, 'linewidth', 1.25); 
% xlabel('Time (Second)'); ylabel('Absolute Echo Density (# Echos / Second)', 'fontsize', fontsize);
% title('Sample RIR', 'fontsize', fontsize + 1); set(gca, 'fontsize', fontsize - 1);
% 
% sh_plt(C, 'horizontal', is_real, 'disp_phase', false, 'dB_lim', [-120, -10], 't', t, 'disp_xaxis_ker_size', 1024);
% sh_plt(C(:, ceil(0.117167 * Fs)), 'mercator', is_real);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate sample absolute echo density curve and directional pdf

% rng(239);
% Fs = 48000;
% duration = 1;
% M = Fs * duration; %Number of samples
% max_odr = 10;
% aed = logspace(log10(10/Fs), log10(2), M) * Fs;
% is_real = true;

% C_pdf = sh_enc_rbf('SqExp', 6, pi/2, 0, 0.5, is_real);
% [C, t] = sh_rand_pp(max_odr, aed, is_real, 'Fs', Fs, 'C_pdf', C_pdf);
% sh_plt(C, 'horizontal', is_real, 'disp_phase', false, 'dB_lim', [-120, 10], 't', t, 'disp_xaxis_ker_size', 1024);

arguments
    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;
    aed (1,:) double {mustBeNonnegative} = 1;
    is_real (1,1) logical = false;

    options.Fs (1,1) double {mustBePositive} = 48000;
    options.sinc_window_halflen (1,1) double {mustBePositive} = 20;
    options.pulse_sample_width (1,1) double {mustBePositive, mustBeInteger} = 1;
    
    options.randomize_phase (1,1) logical = false;
    options.normalize_energy (1,1) logical = true;

    options.direct_unity_first_pulse (1,1) logical = false;

    options.C_pdf (:,:) double = [];

    % Plotting options
    options.enable_disp (1,1) logical = false;
    options.disp_dB_lim (1,2) double = [-120, 0];
    options.disp_xaxis_ker_size (1,1) double {mustBePositive, mustBeInteger} = 1024;
    options.disp_theta_phi (1,2) double = [pi/2, 0];

end

if ~isempty(options.C_pdf)
    assert(sh_pdf_check(options.C_pdf, is_real), 'C_pdf expansion not a density');
end

M = numel(aed);

C = complex(zeros((max_odr + 1)^2, M));

a = options.sinc_window_halflen; %lanczos_kernel half-window size

sigma = 1 ./ sqrt(aed); %Pulse amplitude std. dev over time 

t = 1;
while t <= M
 
    %Draw a sample from exponential distribution   
    tau =  exprnd(options.Fs / aed(floor(t))); %p(tau) = 1/mu * exp(-tau/mu)
    t = t + tau;

    %Sinc kernel
    t_idx = min(max(floor(t-a):ceil(t+a), 1), M); %Write indices
    if t <= M && ~isempty(t_idx)
        
        %Generate pulse
        if options.direct_unity_first_pulse % Force first pulse to be unity
            amp = 1;      
        else
            amp = normrnd(0, 1);
        end
        
        if ~options.randomize_phase
            amp = abs(amp);
        end

        if options.normalize_energy
            amp = amp * sigma(floor(t));
        end

        y_pulse = lanczos_kernel(t_idx - t, a) * amp;

        if options.pulse_sample_width > 1 %Pulse is convolved with non-unity hann window
            y_pulse_conv = circshift(conv(y_pulse, hann(options.pulse_sample_width)'), -floor(options.pulse_sample_width/2));
            y_pulse = y_pulse_conv(1:numel(t_idx));
        end

        %Generate randomized spherical coordinate
        if options.direct_unity_first_pulse %Force theta = pi/2, phi = 0
            theta = pi/2;
            phi = 0;
        else
            if isempty(options.C_pdf)
                theta = acos(2 * rand(1) - 1);
                phi = rand(1) * 2 * pi;
            else
                [theta, phi] = sh_pdf_sample(options.C_pdf(:, min(floor(t), size(options.C_pdf, 2))), 1, "inverseTransform", is_real); 
            end
        end
        
        C_t = sh_enc_proj(max_odr, theta, phi, is_real);
        C_t = C_t / sh_dec(C_t, theta, phi, is_real);  % Normalize Dirac function's peak to unity

        C(:, t_idx) = C(:, t_idx) + C_t * y_pulse;

        options.direct_unity_first_pulse = false; % no longer first-pulse
    end   
end
t = (0:M-1) / options.Fs;

% Plotting
h_fig = [];
if options.enable_disp  && coder.target('MATLAB')

    % Decode RIR in direction of options.disp_theta_phi
    y_theta_phi = real(sh_dec(C,      options.disp_theta_phi(1), options.disp_theta_phi(2), is_real))'; %Decode at options.disp_theta_phi
    y0          = real(sh_dec(C(1,:), options.disp_theta_phi(1), options.disp_theta_phi(2), is_real))'; %Decode 0th order SH

    fontsize = 16;
    h_fig{1} = figure;
    h_fig{1}.Position = [100, 100, 1200, 450];

    tiledlayout(1,2);
    nexttile;
    yyaxis left;  plot(t, y0, 'linewidth', 1.5);
    grid on; axis tight; ylabel('Amplitude', 'fontsize', fontsize);
    yyaxis right;  plot(t, aed, 'linewidth', 1.25); 
    xlabel('Time (Second)'); ylabel('Absolute Echo Density (# Echos / Second)', 'fontsize', fontsize);
    title(['Sample RIR 0th Order Expansion'], ...
        'fontsize', fontsize + 1);
    set(gca, 'fontsize', fontsize - 1);

    nexttile
    yyaxis left;  plot(t, y_theta_phi, 'linewidth', 1.5);
    grid on; axis tight; ylabel('Amplitude', 'fontsize', fontsize);
    yyaxis right;  plot(t, aed, 'linewidth', 1.25); 
    xlabel('Time (Second)'); ylabel('Absolute Echo Density (# Echos / Second)', 'fontsize', fontsize);
    title(['Sample RIR (\theta, \phi) = (', num2str(rad2deg(options.disp_theta_phi(1))), ', ', num2str(rad2deg(options.disp_theta_phi(2))), ')'], ...
        'fontsize', fontsize + 1);
    set(gca, 'fontsize', fontsize - 1);

    % Plot expansions projected onto horizontal plane
    h_tmp = sh_plt(C, 'horizontal', is_real, 'disp_phase', false, 'dB_lim',  options.disp_dB_lim, 't', t, 'disp_xaxis_ker_size', options.disp_xaxis_ker_size);
    h_fig{2} = h_tmp{1};
end