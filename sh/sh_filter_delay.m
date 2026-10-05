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
% C = sh_freqz(sh_enc_proj(3, pi/2, 0), w); % Frequency response of delay-less SH projection
% H = sh_dec(C, pi/2, 0, is_real);
% h = real(ifft(H));
% D = sh_filter_delay(C, w, 100);
% H_del = sh_dec(D, pi/2, 0, is_real);
% h_del = real(ifft(H_del));
% 
% figure; plot(1:N_taps, h, 'ro-', 1:N_taps, h_del, 'b*-'); 
% legend('Original', 'Delayed'); grid on; axis tight;

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
    w (1,:) double = 0;
    d (1,1) double = 0;
end

D = bsxfun(@times, C, exp(-1i * d * w));