function plot_ft_exp_T60_example
%Plot example of exp conv fitted to sample T60 target

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(4123);

Fs = 48000;

T = 1.1;
M = ceil(T * Fs);
h = randn(1, M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify filter g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hz_uniform = linspace(0, Fs/2, 256)';

hz_fc = 4000;
%hz_fc = 8000;
%hz_fc = 12000;

%alpha = 1;
alpha = 1.5;

beta = 1;
%beta = 2;
%beta = 3;

RT60_sec        = mu_lpf(hz_uniform * 2 * pi, alpha, beta, hz_fc * 2 * pi);
RT60_wav_sec    = mu_lpf(hz_uniform * 2 * pi, alpha, 1 * beta, 1 * hz_fc * 2 * pi);
%RT60_wav_sec = mu_pow(hz_uniform * 2 * pi, 0.5 * alpha, -0.1 * beta);

%Plot T60 target
fontsize = 18;
figure; semilogx(hz_uniform, RT60_sec, hz_uniform, RT60_wav_sec, 'LineWidth', 1.5); 
grid on; axis tight; 
xlabel('Hz', 'fontsize', fontsize); ylabel('T60 (Second)', 'fontsize', fontsize);
title('T60 Target', 'fontsize', fontsize + 1);
legend('Noise', 'wav');
set(gca, 'fontsize', fontsize - 1);

%Fit filter g
fit_method = 'constr_min_phase_ls';

[g, g_minphase] = ft_exp_design(4, 'RT60', 'enable_disp', true, 'Fs', Fs, 'RT60_sec', RT60_sec, 'fit_method', fit_method);
g

[g_wav, g_wav_minphase] = ft_exp_design(4, 'RT60', 'enable_disp', true, 'Fs', Fs, 'RT60_sec', RT60_wav_sec, 'fit_method', fit_method);
g_wav

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time-vary convolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_FFT = 4096;
dB_lim = [-60, 0] - 20;
freq_lim = [50, Fs/2];

%position = [100, 100, 800, 400];
position = [100, 100, 800, 360];

[f_exp, f_fig_exp] = ft_exp_conv_opt(h, g, 'enable_disp', true, 'Fs', Fs, 'N_FFT', N_FFT, 'dB_lim', dB_lim, 'freq_lim', freq_lim, 'fontsize', 20, 'disp_position', position);
%[f_rec, f_fig_rec] = ft_rec_conv_opt(h, g, 'enable_disp', true, 'Fs', Fs,  'N_FFT', N_FFT, 'dB_lim', dB_lim, 'freq_lim', freq_lim, 'fontsize', 20, 'disp_position', position);
%err_NMSE(f_exp, f_rec)

%File
[h_wav] = audioread('dat/sample_RIR.wav');
h_wav = h_wav * db2mag(30);
[f_wav_exp, f_wav_fig_exp] = ft_exp_conv_opt(h_wav, g_wav, 'enable_disp', true, 'Fs', Fs, 'N_FFT', N_FFT, 'dB_lim', dB_lim, 'freq_lim', freq_lim, 'fontsize', 20, 'disp_position', position);
xlim(f_wav_fig_exp.Children.Children(2), [-inf, numel(h) / Fs]);


%Export figures
out_dir = 'figs/figs_t60';
if ~isfolder(out_dir)
    mkdir(out_dir);
end

% exportgraphics(f_fig_exp, fullfile(out_dir, 'sample_exp_IR.png'));
% exportgraphics(f_wav_fig_exp, fullfile(out_dir, 'sample_wav_exp_IR.png'));
