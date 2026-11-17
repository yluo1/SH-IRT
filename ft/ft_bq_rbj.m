function [b, a] = ft_bq_rbj(mode, fc, Q, gain_dB, Fs, options)
% Biquadratic cookbook filters from Robert-Bristow-Johnson (RBJ) with 
% additional complementary phase-matched variations

% Phase-matched filter table:
% LPF_PM            ~   HPF
% HPF_PM            ~   LPF
% LS_PM             ~   HS_PM
% PEAK_PM(gain_dB)  ~   PEAK_PM(-gain_dB)
% BPF * BPF         ~   LPF * HPF
% BPF * BPF         ~  -APF

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%mode:       String, filter name {'LPF', 'LPF_PM', 'HPF', 'HPF_PM', 'LS', 'LS_PM', 'HS', 'HS_PM', 'PEAK', 'PEAK_PM', 'BPF', 'NOTCH', 'APF'}

%fc:         Center frequency
%Q:          Bandwidth
%gain_dB:    Filter gain (dB)
%Fs:         Sampling rate

%options:               Struct
%options.enable_disp:   Logical, if true, plot frequency response

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%b:     [1 x 3] Numerator coefficients
%a:     [1 x 3] Denominator coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate and plot phase-matched filters and frequency responses

% fc = 1000;
% Q = sqrt(2)/2;
% gain_dB = 6;
% Fs = 48000;

% % Low-pass
% [b_LPF_PM, a_LPF_PM] = ft_bq_rbj('LPF_PM', fc, Q, gain_dB, Fs);
% [b_HPF,    a_HPF]    = ft_bq_rbj('HPF', fc, Q, gain_dB, Fs);
% fvtool(b_LPF_PM, a_LPF_PM, b_HPF, a_HPF, 'fs', Fs);

% % High-pass
% [b_HPF_PM, a_HPF_PM] = ft_bq_rbj('HPF_PM', fc, Q, gain_dB, Fs);
% [b_LPF,    a_LPF]    = ft_bq_rbj('LPF', fc, Q, gain_dB, Fs);
% fvtool(b_HPF_PM, a_HPF_PM, b_LPF, a_LPF, 'fs', Fs);

% % Shelving
% [b_LS_PM, a_LS_PM] = ft_bq_rbj('LS_PM', fc, Q, gain_dB, Fs);
% [b_HS_PM, a_HS_PM] = ft_bq_rbj('HS_PM', fc, Q, gain_dB, Fs);
% fvtool(b_LS_PM, a_LS_PM, b_HS_PM, a_HS_PM, 'fs', Fs);

% % Peaking
% [b_PEAK_PM, a_PEAK_PM] = ft_bq_rbj('PEAK_PM', fc, Q, gain_dB, Fs);
% [b_CUT_PM, a_CUT_PM]   = ft_bq_rbj('PEAK_PM', fc, Q, -gain_dB, Fs);
% fvtool(b_PEAK_PM, a_PEAK_PM, b_CUT_PM, a_CUT_PM, 'fs', Fs);

% % Band-pass and all-pass
% [b_BPF, a_BPF] = ft_bq_rbj('BPF', fc, Q, gain_dB, Fs);
% [b_APF, a_APF] = ft_bq_rbj('APF', fc, Q, gain_dB, Fs);
% fvtool([b_BPF, a_BPF; b_BPF, a_BPF], [b_LPF, a_LPF; b_HPF, a_HPF], 'fs', Fs);
% fvtool([b_BPF, a_BPF; b_BPF, a_BPF], -b_APF, a_APF, 'fs', Fs);

arguments
    mode (1,:) char {mustBeMember(mode, {'LPF', 'LPF_PM', 'HPF', 'HPF_PM', 'BPF', 'NOTCH', 'APF', 'PEAK', 'PEAK_PM', 'LS', 'LS_PM', 'HS', 'HS_PM'})} = 'LPF';
    fc (1,1) double {mustBePositive} = 1000;
    Q (1,1) double {mustBePositive} = sqrt(2) / 2;
    gain_dB (1,1) double = 0;
    Fs (1,1) double {mustBePositive} = 48000;

    options.enable_disp (1,1) logical = false;
end

A       = 10^(gain_dB/40);
w0      = 2 * pi * fc / Fs;
alpha   = sin(w0) / (2 * Q);

if strcmp(mode, 'LPF')

    b0 =  (1 - cos(w0)) / 2;
    b1 =   1 - cos(w0);
    b2 =  (1 - cos(w0)) / 2;
    a0 =   1 + alpha;
    a1 =  -2 * cos(w0);
    a2 =   1 - alpha;
            
elseif strcmp(mode, 'LPF_PM')

    % -b
    b0 = -(1 - cos(w0)) / 2;
    b1 = -(1 - cos(w0));
    b2 = -(1 - cos(w0)) / 2;
    a0 =   1 + alpha;
    a1 =  -2 * cos(w0);
    a2 =   1 - alpha;
            
elseif strcmp(mode, 'HPF')

    b0 =  (1 + cos(w0)) / 2;
    b1 = -(1 + cos(w0));
    b2 =  (1 + cos(w0)) / 2;
    a0 =   1 + alpha;
    a1 =  -2 * cos(w0);
    a2 =   1 - alpha;
           
elseif strcmp(mode, 'HPF_PM')

    % -b
    b0 = -(1 + cos(w0)) / 2;
    b1 =  (1 + cos(w0));
    b2 = -(1 + cos(w0)) / 2;
    a0 =   1 + alpha;
    a1 =  -2 * cos(w0);
    a2 =   1 - alpha;

elseif strcmp(mode, 'LS')

    b0 =    A * ( (A+1) - (A-1) * cos(w0) + 2 * sqrt(A) * alpha );
    b1 =  2*A * ( (A-1) - (A+1) * cos(w0) );
    b2 =    A * ( (A+1) - (A-1) * cos(w0) - 2 * sqrt(A) * alpha );
    a0 =          (A+1) + (A-1) * cos(w0) + 2 * sqrt(A) * alpha;
    a1 =   -2 * ( (A-1) + (A+1) * cos(w0) );
    a2 =          (A+1) + (A-1) * cos(w0) - 2 * sqrt(A) * alpha;
            
elseif strcmp(mode, 'LS_PM')

    % flip(b)
    b0 =    A * ( (A+1) - (A-1) * cos(w0) - 2 * sqrt(A) * alpha );
    b1 =  2*A * ( (A-1) - (A+1) * cos(w0) );
    b2 =    A * ( (A+1) - (A-1) * cos(w0) + 2 * sqrt(A) * alpha );
    a0 =          (A+1) + (A-1) * cos(w0) + 2 * sqrt(A) * alpha;
    a1 =   -2 * ( (A-1) + (A+1) * cos(w0) );
    a2 =          (A+1) + (A-1) * cos(w0) - 2 * sqrt(A) * alpha;

elseif strcmp(mode, 'HS')
    
    b0 =    A * ( (A+1) + (A-1) * cos(w0) + 2 * sqrt(A) * alpha );
    b1 = -2*A * ( (A-1) + (A+1) * cos(w0) );
    b2 =    A * ( (A+1) + (A-1) * cos(w0) - 2 * sqrt(A) * alpha );
    a0 =          (A+1) - (A-1) * cos(w0) + 2 * sqrt(A) * alpha;
    a1 =    2 * ( (A-1) - (A+1) * cos(w0) );
    a2 =          (A+1) - (A-1) * cos(w0) - 2 * sqrt(A) * alpha;
            
elseif strcmp(mode, 'HS_PM')
    
    % flip(b)
    b0 =    A * ( (A+1) + (A-1) * cos(w0) - 2 * sqrt(A) * alpha );
    b1 = -2*A * ( (A-1) + (A+1) * cos(w0) );
    b2 =    A * ( (A+1) + (A-1) * cos(w0) + 2 * sqrt(A) * alpha );
    a0 =          (A+1) - (A-1) * cos(w0) + 2 * sqrt(A) * alpha;
    a1 =    2 * ( (A-1) - (A+1) * cos(w0) );
    a2 =          (A+1) - (A-1) * cos(w0) - 2 * sqrt(A) * alpha;

elseif strcmp(mode, 'PEAK')

    b0 =   1 + alpha * A;
    b1 =  -2 * cos(w0);
    b2 =   1 - alpha * A;
    a0 =   1 + alpha / A;
    a1 =  -2 * cos(w0);
    a2 =   1 - alpha / A;

elseif strcmp(mode, 'PEAK_PM')

    % flip(b)
    b0 =   1 - alpha * A;
    b1 =  -2 * cos(w0);
    b2 =   1 + alpha * A;
    a0 =   1 + alpha / A;
    a1 =  -2 * cos(w0);
    a2 =   1 - alpha / A;

elseif strcmp(mode, 'BPF')

    b0 =   alpha;
    b1 =   0;
    b2 =  -alpha;
    a0 =   1 + alpha;
    a1 =  -2 * cos(w0);
    a2 =   1 - alpha;

elseif strcmp(mode, 'NOTCH')

    b0 =   1;
    b1 =  -2 * cos(w0);
    b2 =   1;
    a0 =   1 + alpha;
    a1 =  -2 * cos(w0);
    a2 =   1 - alpha;

elseif strcmp(mode, 'APF')

    b0 =   1 - alpha;
    b1 =  -2 * cos(w0);
    b2 =   1 + alpha;
    a0 =   1 + alpha;
    a1 =  -2 * cos(w0);
    a2 =   1 - alpha;

else
    
    error('unknown filter type');
end

% Normalize by a0
b = [b0, b1, b2] / a0;
a = [a0, a1, a2] / a0;

% Plotting
if options.enable_disp && coder.target("MATLAB")
    fvtool(b, a, 'Fs', Fs);
end