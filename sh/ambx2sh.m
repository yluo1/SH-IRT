function D = ambx2sh(C)
%Convert AmbiX to N3D spherical harmonics expansion

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:            [(P + 1)^2 x M] AmbiX coefficients (max order P of M number of functions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%D:            [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Read AmbiX file and decode at cardinal directions

% [C_amb, Fs] = audioread('dat/voice_noise_ambx.wav');
% C_amb = C_amb';
% t = (0:(size(C_amb, 2)-1)) / Fs;

% C = ambx2sh(C_amb(:, :));
% f_left = real(sh_dec(C, pi/2, deg2rad(-90), false, 'eval'));
% f_right = real(sh_dec(C, pi/2, deg2rad(90), false, 'eval'));
% f_top = real(sh_dec(C, 0, 0, false, 'eval'));
% f_bot = real(sh_dec(C, pi, 0, false, 'eval'));
% f_center = real(sh_dec(C, pi/2, 0, false, 'eval'));
% f_back = real(sh_dec(C, pi/2, pi, false, 'eval'));

% fontsize = 16;
% figure; plot(t, f_center, t, f_back, ':', 'linewidth', 1.5); grid on; axis tight;
% h_lg = legend('Center', 'Back', 'location', 'best'); set(h_lg, 'fontsize', fontsize - 1);
% xlabel('Time (Second)', 'fontsize', fontsize);
% ylabel('Amplitude', 'fontsize', fontsize);
% title('AmbiX Directional Decoding', 'fontsize', fontsize + 1);
% set(gca, 'fontsize', fontsize - 1);

% sh_plt(C(:, 24000), 'mercator', false);


arguments
    C (:,:) double = 0;
end

[P, M] = size(C);
P = sqrt(P) - 1;
assert(P - floor(P) == 0, 'Invalid size C');

D_real = C;
for l = 0:P
    idx = (l^2 + 1) : (l+1)^2;
    D_real(idx, :) = C(idx, :) .* sqrt(2 * l + 1);              
end

D = sh_re2cpx(D_real);
