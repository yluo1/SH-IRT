function plot_ft_exp_design_example
%Plot exponetiated filter design examples for target RT60 at varying number of filter taps

%Author: Yuancheng Luo, 2026

Fs = 48000;

RT60_sec = [4 0.25 4 4 0.5]; %Target RT60 from DC to Nyquist with overshoot 
%RT60_sec = [4 0.5 2 2 0.5] * 2;
%RT60_sec = [5 1.5 4 4 2.5 1 0.3 0.25 0.24 0.2 0.22];

tol0 = 1e-6;

%g_RT60 = ft_exp_design(9, 'RT60',  'enable_disp', true, 'Fs', Fs, 'RT60_sec', RT60_sec, 'tol0', tol0);
%g_RT60 = ft_exp_design(5, 'RT60',  'enable_disp', true, 'Fs', Fs, 'RT60_sec', RT60_sec, 'tol0', tol0);

H_tgt_dB_DC_NQ  = -60 ./ (Fs * RT60_sec); %Target magnitude dB
H_tgt_abs_DC_NQ = db2mag(H_tgt_dB_DC_NQ);                 %Target magnitude modulus  

X_mag = [H_tgt_abs_DC_NQ, fliplr(H_tgt_abs_DC_NQ(2:end-1))];

[g_tap_8, g_minphase, err, lambda]  = ft_bnd_minphase(X_mag(:), 8, 'ub', 1 - tol0, 'enable_disp', false);
[g_tap_6]                           = ft_bnd_minphase(X_mag(:), 6, 'ub', 1 - tol0, 'enable_disp', false);
[g_tap_4]                           = ft_bnd_minphase(X_mag(:), 4, 'ub', 1 - tol0, 'enable_disp', false);
[g_tap_3]                           = ft_bnd_minphase(X_mag(:), 3, 'ub', 1 - tol0, 'enable_disp', false);
[g_tap_2]                           = ft_bnd_minphase(X_mag(:), 2, 'ub', 1 - tol0, 'enable_disp', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RT60_freq = linspace(0, Fs/2, numel(RT60_sec));
% RT60_freq(1) = 20; %For display purposes
% disp_freq = logspace(log10(20), log10(Fs/2), 512);

RT60_freq = linspace(0, Fs/2, numel(RT60_sec));
disp_freq = linspace(0, Fs/2, 512);


H_minphase = freqz(g_minphase, 1, disp_freq, Fs);
H_g_tap_8 = freqz(g_tap_8, 1, disp_freq, Fs);
H_g_tap_6 = freqz(g_tap_6, 1, disp_freq, Fs);
H_g_tap_4 = freqz(g_tap_4, 1, disp_freq, Fs);
H_g_tap_3 = freqz(g_tap_3, 1, disp_freq, Fs);
H_g_tap_2 = freqz(g_tap_2, 1, disp_freq, Fs);

grpdel_minphase = grpdelay(g_minphase, 1, disp_freq, Fs);
grpdel_g_tap_8 = grpdelay(g_tap_8, 1, disp_freq, Fs);
grpdel_g_tap_6 = grpdelay(g_tap_6, 1, disp_freq, Fs);
grpdel_g_tap_4 = grpdelay(g_tap_4, 1, disp_freq, Fs);
grpdel_g_tap_3 = grpdelay(g_tap_3, 1, disp_freq, Fs);
grpdel_g_tap_2 = grpdelay(g_tap_2, 1, disp_freq, Fs);
grpdel_tgt     = grpdelay(g_minphase, 1, RT60_freq, Fs);

%0.25 sec exponentiation
sample_exp = ceil(Fs * 0.25);

H_minphase = H_minphase.^sample_exp;
H_g_tap_8 = H_g_tap_8.^sample_exp;
H_g_tap_6 = H_g_tap_6.^sample_exp;
H_g_tap_4 = H_g_tap_4.^sample_exp;
H_g_tap_3 = H_g_tap_3.^sample_exp;
H_g_tap_2 = H_g_tap_2.^sample_exp;
H_tgt_abs_DC_NQ = H_tgt_abs_DC_NQ.^sample_exp;

grpdel_minphase = grpdel_minphase * sample_exp;
grpdel_g_tap_8 = grpdel_g_tap_8 * sample_exp;
grpdel_g_tap_6 = grpdel_g_tap_6 * sample_exp;
grpdel_g_tap_4 = grpdel_g_tap_4 * sample_exp;
grpdel_g_tap_3 = grpdel_g_tap_3 * sample_exp;
grpdel_g_tap_2 = grpdel_g_tap_2 * sample_exp;
grpdel_tgt = grpdel_tgt * sample_exp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fontsize = 16;
backgroundAlpha = 1;
markersize = 12;

h_f = figure;
h_f.Position = [100, 100, [560, 420] * 1];

tiledlayout(2, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Magnitude plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexttile;
plot(   RT60_freq, mag2db(abs(H_tgt_abs_DC_NQ)), 'r*', ...
            disp_freq, mag2db(abs(H_minphase)), ...
            disp_freq, mag2db(abs(H_g_tap_8)), ...
            disp_freq, mag2db(abs(H_g_tap_6)), ...
            disp_freq, mag2db(abs(H_g_tap_4)), ...
            disp_freq, mag2db(abs(H_g_tap_2)), ...
            'linewidth', 1.5, 'MarkerSize', markersize); 

grid on; axis tight;

xlabel('Frequency (Hz)', 'fontsize', fontsize);
ylabel('Magnitude (dB)', 'fontsize', fontsize);
title('Filter Frequency Response: 0.25 Second Exponentiation', 'fontsize', fontsize + 1);

h_lg = legend({ 'Target Responses', ...
                'Minimum Phase', ...
                'Fit N = 8', ...
                'Fit N = 6', ...
                'Fit N = 4', ...
                'Fit N = 2', ...
             }, 'location', 'southeast', 'NumColumns', 2, 'BackgroundAlpha', backgroundAlpha);
set(h_lg, 'fontsize', fontsize - 2);
set(gca, 'fontsize', fontsize - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Phase plot
nexttile;

plot(   RT60_freq, grpdel_tgt, 'r*', ...
            disp_freq, grpdel_minphase, ...
            disp_freq, grpdel_g_tap_8, ...
            disp_freq, grpdel_g_tap_6, ...
            disp_freq, grpdel_g_tap_4, ...
            disp_freq, grpdel_g_tap_2, ...
            'linewidth', 1.5, 'MarkerSize', markersize); 

grid on; axis tight;
xlabel('Frequency (Hz)', 'fontsize', fontsize);
ylabel('Group Delay (Sample)', 'fontsize', fontsize);

h_lg = legend({ 'Target Responses', ...
                'Minimum Phase', ...
                'Fit N = 8', ...
                'Fit N = 6', ...
                'Fit N = 4', ...
                'Fit N = 2', ...
             }, 'location', 'southeast', 'NumColumns', 2, 'BackgroundAlpha', backgroundAlpha);
set(h_lg, 'fontsize', fontsize - 2);
set(gca, 'fontsize', fontsize - 1);

%Export figures
out_dir = 'figs';
if ~isfolder(out_dir)
    mkdir(out_dir);
end
exportgraphics(h_f, fullfile(out_dir, 'sample_filter_fit.png'));