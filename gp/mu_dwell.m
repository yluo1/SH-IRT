function mu = mu_dwell(omega, mode, vol)
% Mean T60 times as a function of room type, cubic volume, and frequency.
% Paper reference: Díaz, C., & Pedrero, A. (2007). The reverberation time and 
% equivalent sound absorption area of rooms in dwellings. Noise & Vibration Worldwide, 38(6), 12-20.

% T60 augmented with linear extrapolated frequency extension of +-2 octave  outside [125, 4000] Hz
% Piecewise spline interpolation and extrapolation for specified frequencies

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%omega:     [N x M] Angular frequency (freq * 2 * pi)
%mode:      String, dwelling type {'furnished', 'unfurnished'}
%vol:       Scalar, cubic volume (vol m^3)
%                   10 <= vol <= 100  for mode = 'furnished'
%                   10 <= vol <= 60   for mode = 'unfurnished'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%mu:        [N x M] T60 (seconds) means

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Plot sample T60 for furnished rooms

% N = 256;
% freq = logspace(log10(10), log10(24000), N)';
% omega = 2 * pi * freq;
% mode = 'furnished';
% mu_T60_V12 = mu_dwell(omega, mode, 12);
% mu_T60_V25 = mu_dwell(omega, mode, 25);
% mu_T60_V50 = mu_dwell(omega, mode, 50);
% mu_T60_V75 = mu_dwell(omega, mode, 75);
% mu_T60_V100 = mu_dwell(omega, mode, 100);

% fontsize = 16;
% figure;
% semilogx(freq, [mu_T60_V12, mu_T60_V25, mu_T60_V50, mu_T60_V75, mu_T60_V100], 'linewidth', 1.5); 
% h_lg = legend('Vol = 12 m^3', 'Vol = 25 m^3', 'Vol = 50 m^3', 'Vol = 75 m^3', 'Vol = 100 m^3', 'location', 'best'); set(h_lg, 'fontsize', fontsize - 1);
% grid on; axis tight; set(gca, 'fontsize', fontsize - 1);
% xlabel('Frequency (Hz)', 'fontsize', fontsize);
% ylabel('T60 (Seconds)', 'fontsize', fontsize);
% title(['T60 Dwelling ', mode], 'fontsize', fontsize + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Plot sample T60 for unfurnished rooms

% N = 256;
% freq = logspace(log10(10), log10(24000), N)';
% omega = 2 * pi * freq;
% mode = 'unfurnished';
% mu_T60_V12 = mu_dwell(omega, mode, 12);
% mu_T60_V25 = mu_dwell(omega, mode, 25);
% mu_T60_V50 = mu_dwell(omega, mode, 50);
% 
% fontsize = 16;
% figure;
% semilogx(freq, [mu_T60_V12, mu_T60_V25, mu_T60_V50], 'linewidth', 1.5); 
% h_lg = legend('Vol = 12 m^3', 'Vol = 25 m^3', 'Vol = 50 m^3',  'location', 'best'); set(h_lg, 'fontsize', fontsize - 1);
% grid on; axis tight; set(gca, 'fontsize', fontsize - 1);
% xlabel('Frequency (Hz)', 'fontsize', fontsize);
% ylabel('T60 (Seconds)', 'fontsize', fontsize);
% title(['T60 Dwelling ', mode], 'fontsize', fontsize + 1);

arguments
    omega (:,1) double {mustBeNonnegative} = 0;
    mode   (1,:) char {mustBeMember(mode, {'furnished', 'unfurnished'})} = 'furnished';
    vol (1,1) double {mustBePositive} = 25;
end

freq_in = omega / (2*pi);

if strcmp(mode, 'furnished')

    freq = [125, 250, 500, 1000, 2000, 4000];
    T_cell = { ...
        @(V)  7.10E-6 * V^2 + 0.0026 * V + 0.452, ...
        @(V) -6.10E-6 * V^2 + 0.004  * V + 0.3236, ...
        @(V) -2.10E-5 * V^2 + 0.0048 * V + 0.255, ...
        @(V) -2.10E-5 * V^2 + 0.0046 * V + 0.2239, ...
        @(V) -2.10E-5 * V^2 + 0.0048 * V + 0.2013, ...
        @(V) -2.10E-5 * V^2 + 0.0042 * V + 0.1929, ...
        };

    assert(all(vol >= 10 & vol <= 100), 'vol must be between [10, 100]')

elseif strcmp(mode, 'unfurnished')

    freq = [125, 250, 500, 1000, 2000, 4000];
    T_cell = { ...
        @(V) -9.10E-4 * V^2 + 0.0875 * V + 0.0824, ...
        @(V) -9.10E-4 * V^2 + 0.0863 * V - 0.0309, ...
        @(V) -6.10E-4 * V^2 + 0.0624 * V + 0.1211, ...
        @(V) 0.0241 * V + 0.3949, ...
        @(V) 0.0219 * V + 0.2977, ...
        @(V) 0.0171 * V + 0.2934, ...
        };

    assert(all(vol >= 10 & vol <= 60), 'vol must be between [10, 60]')
    
else
    error('Unknown mode');
end

% Augmented frequencies with +-2 octave extension with linear extrapolation
freq_aug = [freq(1)/4, freq(1)/2, freq, freq(end)*2, freq(end)*4]; 
log_mu_aug = interp1( log(1 + freq), log(abs(cellfun(@(f) f(vol), T_cell)) ), log(1 + freq_aug), 'linear', 'extrap');

% Evaluate at arbitrary frequency. Extrapolate with piecewise spline derivatives beyond augmented frequencies 
log_mu = interp1( log(1 + freq_aug), log_mu_aug, log(freq_in + 1), 'pchip', 'extrap');
mu = exp(log_mu);
