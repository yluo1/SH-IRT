function D = sh_filter(C, varargin)
%Filter spherical harmonic expansion

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:         [(P + 1)^2 x M] SH coefficients

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
%Sample usage:  Filter SH random field

% Fs    = 48000;
% M_w   = 256;
% freq  = logspace(log10(100), log10(Fs / 2), M_w);
% w     = freq / Fs * 2 * pi;

% fc = 2000;
% [b,a]     = butter(4,fc/(Fs/2));
% [z,p,k]   = butter(4,fc/(Fs/2));
% sos       = zp2sos(z,p,k);

% rng(234);
% max_odr = 3;
% num_func = 1024;

% is_real = true;
% C = sh_rand(max_odr, num_func, is_real, true);
% C_butter_ba  =  sh_filter(C, b, a);
% C_butter_sos =  sh_filter(C, sos);

% D     = sh_freqz(C, w);
% D_ba  = sh_freqz(C_butter_ba, w);
% D_sos = sh_freqz(C_butter_sos, w);

% dB_lim = [-10, 70];
% sh_plt(D, 'horizontal', is_real, 'hz', freq, 'title_name', 'Real Random Field', 'dB_lim', dB_lim);
% sh_plt(D_ba, 'horizontal', is_real, 'hz', freq, 'title_name', 'Real Random Field * LPF (b/a)', 'dB_lim', dB_lim);
% sh_plt(D_sos, 'horizontal', is_real, 'hz', freq, 'title_name', 'Real Random Field * LPF SOS', 'dB_lim', dB_lim);

% err = norm(D_ba - D_sos)

arguments
    C (:,:) double = 0;
end

arguments (Repeating)
    varargin
end

N_varagin = numel(varargin);

if N_varagin == 0       %Identity

    D = C;

elseif N_varagin <= 1   %Second-order sections (SOS) matrix
    
    sos = varargin{1};
    D = sosfilt(sos, C, 2);

elseif N_varagin <= 2    %Numerator / Denominator filter cofficients b, a

    b = varargin{1}(:).';
    a = varargin{2}(:).';    
    D = filter(b, a, C, [],2); 
    
elseif N_varagin <= 3    %Mixed coefficient vectors and matrix b, a, SOS
    
    b = varargin{1}(:).';
    a = varargin{2}(:).';
    sos = varargin{3};
    D = sosfilt(sos, filter(b, a, C, [], 2), 2); 
       
end