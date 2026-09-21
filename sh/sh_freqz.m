function D = sh_freqz(C, w)
%Frequency response of spherical harmonic expansion

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:         [(P + 1)^2 x M] SH coefficients
%w:         [1 x M_w] Angular frequency (radians / sample)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%D:          [(P + 1)^2 x M] SH coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Get frequency response of SH random field

% Fs    = 48000;
% M_w   = 256;
% freq  = logspace(log10(100), log10(Fs / 2), M_w);
% w     = freq / Fs * 2 * pi;

% rng(234);
% max_odr = 3;
% num_func = 1024;

% is_real = true;
% C = sh_rand(max_odr, num_func, is_real, true);
% C_freqz  =  sh_freqz(C, w);

% dB_lim = [-10, 70];
% sh_plt(C_freqz, 'horizontal', is_real, 'hz', freq, 'title_name', 'Real Random Field', 'dB_lim', dB_lim);

arguments
    C (:,:) double = 0;
    w (1,:) double = 0;
end

[N_C, M] = size(C);
M_w = numel(w);

w_aug = [0; w(:)]';     %Include DC component for freqz

D = complex(zeros(N_C, M_w));
for n = 1:N_C
    D_aug_n = freqz(C(n, :), 1, w_aug);
    D(n, :) = D_aug_n(2:end);
end