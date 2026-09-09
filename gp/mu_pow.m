function mu = mu_pow(omega, alpha, beta)
%Evaluate power-law function of angular frequency:
%mu(omega) = alpha * max(1, omega / (2 * pi)).^(-beta)

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%omega:     [N x M] Angular frequency
%alpha:     Scalar, scaling coefficient, non-negative
%beta:      Scalar, power

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%mu:        [N x M] Means

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Plot sample functions

% N = 256;
% freq = logspace(log10(20), log10(24000), N)';
% omega = 2 * pi * freq;
% alpha = 1;
% beta = 0.25;
% mu = mu_pow(omega, alpha, beta);
% plot_cpx_resp(mu, freq, 'title_name', 'mu pow');

arguments
    omega (:,1) double {mustBeNonnegative} = pi/2;

    alpha   (1,1) double {mustBeNonnegative} = 1;
    beta    (1,1) double = 0.25;
end

mu = alpha * max(1, omega / (2 * pi)).^(-beta);