function [h] = ft_two_tap_FIR(dB_DC, dB_NQ, enable_disp)
%Two-tap min-phase FIR with DC and Nyquist gain specifications

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%dB_DC:         dB gain at DC
%dB_NQ:         dB gain at Nyquist
%enable_disp:   Logical, if true, plot frequency response assuming 48 kHz Fs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h:             [1 x 2] FIR coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage

%h = ft_two_tap_FIR(0, -6, true);

arguments
    dB_DC (1,1) double = 0;
    dB_NQ (1,1) double = 0;
    enable_disp (1,1) logical = false;
end

mag_DC = db2mag(dB_DC);
mag_NQ = db2mag(dB_NQ);

mag_sum  = (mag_DC + mag_NQ) / 2;
mag_diff = (mag_DC - mag_NQ) / 2;

if abs(mag_sum) > abs(mag_diff)
    h = [mag_sum, mag_diff];
else
    h = [mag_diff, mag_sum];
end

%Plotting
if enable_disp
    fvtool(h, 1, 'fs', 48000);
end