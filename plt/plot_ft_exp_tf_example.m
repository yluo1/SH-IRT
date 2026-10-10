function plot_ft_exp_tf_example
% Generate sample transfer function plot of exponentiated self-convolutions g,
% Gaussian white-noise h, and exponentiated filtered f

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate exponentiated filter g
Fs = 8000; %Sample rate
RT60_sec = [0.45, 0.55, 0.5, 0.45, 0.4, 0.36, 0.3, 0.275]; % Target RT60 from DC to Nyquist with overshoot 

tol0 = 1e-6;

[g] = ft_exp_design(numel(RT60_sec), 'RT60',  'enable_disp', false, 'Fs', Fs, 'RT60_sec', RT60_sec, 'tol0', tol0);

hz = logspace(log10(50), log10(Fs/2), 512)';
G = freqz(g, 1, hz, Fs);

% Generate RIR h
rng(451);
h = randn(1, ceil(0.6 * Fs));

% Exponentiated convolution
[f, h_fig] = ft_exp_conv_opt(h, g, 'enable_disp', true, 'Fs', Fs, 'N_FFT', 512);

%Plot g
h_g = plot_cpx_resp([G, G.^(10/1000 * Fs), G.^(100/1000 * Fs)], hz, ...
    'fontsize', 30, 'linewidth', 2.5, 'title_name', 'Filters g(m) Frequency Responses', ...
    'resp_names', {'g(1) @ 0 ms', 'g(80) @ 10 ms', 'g(800) @ 100 ms'});

% Plot f
fontsize = 18;
h_h = figure;
plot((0:numel(h)-1)/ Fs * 1000, h, 'linewidth', 1.25)
xlabel('Time (ms)', 'fontsize', fontsize);
ylabel('Amplitude', 'fontsize', fontsize);
title('Filter h (Gaussian White Noise)', 'fontsize', fontsize + 1);
grid on; axis tight; set(gca, 'fontsize', fontsize - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_dir = 'figs/figs_t60';
if ~isfolder(out_dir)
    mkdir(out_dir);
end
exportgraphics(h_h, fullfile(out_dir, 'exp_h_tf_example.png'));
exportgraphics(h_g, fullfile(out_dir, 'exp_g_tf_example.png'));
exportgraphics(h_fig, fullfile(out_dir, 'exp_f_tf_example.png'));
