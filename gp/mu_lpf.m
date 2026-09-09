function mu = mu_lpf(omega, alpha, beta, omega_fc)
%Evaluate low-pass filter magnitude function of frequency:
%mu(omega) = alpha ./ (1 + (omega ./ omega_fc).^beta);

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%omega:     [N x M] Angular frequency
%alpha:     Scalar, scaling coefficient, non-negative
%beta:      Scalar, power, non-negative
%omega_fc:  Scalar, corner angular frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%mu:        [N x M] Means

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Plot sample functions

% N = 256;
% freq = logspace(log10(20), log10(24000), N)';
% omega = 2 * pi * freq;
% alpha = 1;
% beta = 1;
% freq_fc = [4000, 8000, 16000];
% omega_fc = 2 * pi * freq_fc;
% mu = [mu_lpf(omega, alpha, beta, omega_fc(1)), mu_lpf(omega, alpha, beta, omega_fc(2)), mu_lpf(omega, alpha, beta, omega_fc(3))];
% plot_cpx_resp(mu, freq, 'title_name', 'mu low-pass-filter', 'resp_name', cellstr(string(freq_fc)) );

arguments
    omega    (:,1) double {mustBeNonnegative} = pi/2;

    alpha    (1,1) double {mustBeNonnegative} = 1;
    beta     (1,1) double {mustBeNonnegative} = 0.25;

    omega_fc (1,1) double {mustBePositive} = pi/2;
end

mu = alpha ./ (1 + (omega ./ omega_fc).^beta);