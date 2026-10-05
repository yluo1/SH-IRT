function D = sh_filter_delay(C, w, d)
%Add sample delay d to SH expansion frequency response C at angular frequency w

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:         [(P + 1)^2 x M] SH coefficients at angular frequencies w
%w:         [1 x M] Angular frequency (radians / sample) [0 to 2 * pi]
%d:         Scalar, sample delay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%D:         [(P + 1)^2 x M] SH coefficients delayed by d samples at angular frequencies w

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate SH projection and add delay

% is_real = false;
% N_taps = 128;
% w = linspace(0, 2 * pi, N_taps + 1); w = w(1:end-1);
% t = (0:(N_taps-1))/Fs;
% Fs = 8000;
% omega = w * Fs;
% freq = omega / (2 * pi);

% C = sh_freqz(sh_enc_proj(3, pi/2, 0), w); % Frequency response of delay-less SH projection
% D = sh_filter_delay(C, w, 10/1000 * Fs); % 10 ms

% sh_plt(C, 'cardinal', is_real, 'hz', freq, 'title_name', 'Original');
% sh_plt(D, 'cardinal', is_real, 'hz', freq, 'title_name', 'Delayed');

% C_td = ifft(C, [], 2, 'symmetric');
% sh_plt(C_td, 'cardinal', is_real, 't', t, 'title_name', 'Original');
% D_td = ifft(D, [], 2, 'symmetric');
% sh_plt(D_td, 'cardinal', is_real, 't', t, 'title_name', 'Delayed');

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
    w (1,:) double = 0;
    d (1,1) double = 0;
end

D = bsxfun(@times, C, exp(-1i * d * w));