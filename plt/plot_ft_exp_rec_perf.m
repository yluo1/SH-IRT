function plot_ft_exp_rec_perf
%Evaluate and plot exponential and recursive convolution runtime performance

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(3451);
Fs = 48000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vary FIR h size M
%M_list_sec = [0.125, 0.25, 0.5, 1.0];
M_list_sec = [0.125, 0.25, 0.5, 1.0, 2.0];

M_list = ceil(M_list_sec * Fs);
M = numel(M_list);
h_list = cell(1, M);
for m = 1:M
    h_list{m} = randn(1, M_list(m));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vary FIR g size N
%N_list = [2, 4, 8, 16, 32];
N_list = [2, 4, 8, 16];

N = numel(N_list);
g_list = cell(1, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify filter g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hz_uniform = linspace(0, Fs/2, 256)';

hz_fc = 4000;
%hz_fc = 8000;
%hz_fc = 12000;

alpha = 2;

beta = 1;
%beta = 2;
%beta = 3;

RT60_sec = mu_lpf(hz_uniform * 2 * pi, alpha, beta, hz_fc * 2 * pi);

%Plot T60 target
fontsize = 18;
figure; semilogx(hz_uniform, RT60_sec, 'LineWidth', 1.5); 
grid on; axis tight; 
xlabel('Hz', 'fontsize', fontsize); ylabel('T60 (Second)', 'fontsize', fontsize);
title('T60 Target', 'fontsize', fontsize + 1);
set(gca, 'fontsize', fontsize - 1);


fit_method = 'constr_min_phase_ls';
%fit_method = 'constr_min_phase_minimax';
%fit_method = 'two_tap';

for n = 1:N
    g_list{n} = ft_exp_design(N_list(n), 'RT60', 'enable_disp', false, 'Fs', Fs, 'RT60_sec', RT60_sec, 'fit_method', fit_method);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Runtime performance
max_iter = 20;

exp_conv_time_sec = zeros([M, N]);
rec_conv_time_sec = zeros([M, N]);

f_exp_list = cell([M, N]);
f_rec_list = cell([M, N]);

for m = 1:M
    for n = 1:N

        %ft_exp_conv_opt
        toc_mn = 0;
        for iter = 1:max_iter
            tic;
            f_exp_list{m, n} = ft_exp_conv_opt(h_list{m}, g_list{n});
            toc_mn = toc_mn + toc;
        end
        exp_conv_time_sec(m, n) = toc_mn / max_iter

        %ft_rec_conv_opt
        toc_mn = 0;
        for iter = 1:max_iter
            tic;
            f_rec_list{m, n} = ft_rec_conv_opt(h_list{m}, g_list{n});
            toc_mn = toc_mn + toc;
        end
        rec_conv_time_sec(m, n) = toc_mn / max_iter

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot IR of longest h, g

% N_FFT = 4096;
% dB_lim = [-80, -20];
% ft_exp_conv_opt(h_list{end}, g_list{end}, 'enable_disp', true, 'N_FFT', N_FFT, 'dB_lim', dB_lim);
% ft_rec_conv_opt(h_list{end}, g_list{end}, 'enable_disp', true, 'N_FFT', N_FFT, 'dB_lim', dB_lim);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Runtime plots

fontsize = 18;
linewidth = 1.5;
markersize = 12;
h_runtime = figure;
h_runtime.Position = [100, 100, 560, 380];

loglog(M_list_sec, exp_conv_time_sec, 'o-', 'linewidth', linewidth, 'markersize', markersize); hold on;
set(gca, 'ColorOrderIndex', 1);
loglog(M_list_sec, rec_conv_time_sec, '*--', 'linewidth', linewidth, 'markersize', markersize); 
xlabel('FIR h Size M (Second)', 'fontsize', fontsize);
ylabel('Runtime  (Second)', 'fontsize', fontsize);
title('Runtimes for Exponentiating FIR Size N', 'fontsize', fontsize + 1);
grid on; axis tight;
set(gca, 'fontsize', fontsize - 1);
ylim([-inf, 10]);
xticks(M_list_sec);

legend_str = {};
for n = 1:N
    legend_str{end+1} = ['Exp. N = ', num2str(N_list(n))];
end
for n = 1:N
    legend_str{end+1} = ['Rec. N = ', num2str(N_list(n))];
end
h_lg = legend(legend_str, 'location', 'northwest', 'NumColumns', 2); 
set(h_lg, 'fontsize', fontsize - 2);


%Export figures
out_dir = 'figs';
if ~isfolder(out_dir)
    mkdir(out_dir);
end
exportgraphics(h_runtime, fullfile(out_dir, 'sample_exp_rec_runtime.png'));