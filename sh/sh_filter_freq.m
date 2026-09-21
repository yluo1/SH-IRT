function D = sh_filter_freq(C, w, varargin)
%Frequency-domain filtering of spherical harmonic expansion coefficients

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:         [(P + 1)^2 x M] SH coefficients
%w:         [1 x M] Angular frequency (radians / sample)

%varargin:  Filter coefficients

%Format:        H(SOS)
%varargin{1}    = SOS  [N_SOS x 6] Second-order sections

%Format:        H(b) / H(a):
%varargin{1}    = b    [N_b x 1] or [1 x N_b] numerator coefficients
%varargin{2}    = a    [N_a x 1] or [1 x N_a] denominator coefficients

%Format:        H(SOS_n) * H(b_n) / H(a_n):
%varargin{1}    = b    [N_b x 1] or [1 x N_b] numerator coefficients
%varargin{2}    = a    [N_a x 1] or [1 x N_a] denominator coefficients
%varargin{3}    = SOS  [N_SOS x 6] Second-order sections

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%D:          [(P + 1)^2 x M] SH coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Filter SH piston frequency response

% Fs    = 48000;
% freq  = logspace(log10(100), log10(Fs / 2), 256);
% w     = freq / Fs * 2 * pi;

% fc = 2000;
% [b,a]     = butter(4,fc/(Fs/2));
% [z,p,k]   = butter(4,fc/(Fs/2));
% sos       = zp2sos(z,p,k);

% C = sh_enc_pist_sphere(40, pi/2, 0, 0.25, deg2rad(20), freq, 4);
% C_butter_ba  =  sh_filter_freq(C, w, b, a);
% C_butter_sos =  sh_filter_freq(C, w, sos);

% dB_lim = [-10, 70];
% sh_plt(C, 'horizontal', false, 'hz', freq, 'title_name', 'Piston on Spherical Baffle', 'dB_lim', dB_lim);
% sh_plt(C_butter_ba, 'horizontal', false, 'hz', freq, 'title_name', 'Piston on Spherical Baffle * LPF (b/a)', 'dB_lim', dB_lim);
% sh_plt(C_butter_sos, 'horizontal', false, 'hz', freq, 'title_name', 'Piston on Spherical Baffle * LPF SOS', 'dB_lim', dB_lim);

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
    w (1,:) double = 0;
end

arguments (Repeating)
    varargin
end

N_varagin = numel(varargin);
M = size(C, 2);

if M ~= numel(w)
    error('Size mismatch C, w');
end

w_aug = [0; w(:)]';     %Include DC component for freqz

if N_varagin == 0       %Identity

    D = C;

elseif N_varagin <= 1   %Second-order sections (SOS) matrix
    
    sos = varargin{1};

    if size(sos, 1) == 1
        H = freqz(sos(1:3), sos(4:6), w_aug);
    else
        H = freqz(sos, w_aug);
    end
    D = bsxfun(@times, C, H(2:end));

elseif N_varagin <= 2    %Numerator / Denominator filter cofficients b, a

    b = varargin{1}(:).';
    a = varargin{2}(:).';
    
    H = freqz(b, a, w_aug);
    D = bsxfun(@times, C, H(2:end));
    
elseif N_varagin <= 3    %Mixed coefficient vectors and matrix b, a, SOS
    
    b = varargin{1}(:).';
    a = varargin{2}(:).';
    sos = varargin{3};
        
    if size(SOS, 1) == 1
        H = freqz(sos(1:3), sos(4:6), w_aug);
    else
        H = freqz(sos, w_aug);
    end
    H = H .* freqz(b, a, w_aug);
    D = bsxfun(@times, C, H(2:end));

end


