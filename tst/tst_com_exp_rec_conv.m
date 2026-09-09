function tst_com_exp_rec_conv
%Test commutative properties of exponentiating and recursive convolutions

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(4123);
Fs = 48000;

%T = 1;
T = 0.1;
%T = 0.01;

M = ceil(T * Fs);
h = randn(1, M);

%Compute Exponentiating filter
g_a = ft_two_tap_FIR(-1e-3, -3, false); %Larger differences
g_b = ft_two_tap_FIR(-5e-2, -4, false); %Larger differences

% g_a = ft_exp_design(6, 'RT60',  'enable_disp', false, 'Fs', 48000, 'RT60_sec', [2 0.35 0.4 0.2 0.3]);
% g_b = ft_exp_design(6, 'RT60',  'enable_disp', false, 'Fs', 48000, 'RT60_sec', [1 0.15 0.5 0.4 0.2]);

%f(f(h, g_a), g_b)
f_exp_h_a_b = ft_exp_conv_opt(ft_exp_conv_opt(h, g_a), g_b, 'enable_disp', false);
f_rec_h_a_b = ft_rec_conv_opt(ft_rec_conv_opt(h, g_a), g_b, 'enable_disp', false);

%f(f(h, g_b), g_a)
f_exp_h_b_a = ft_exp_conv_opt(ft_exp_conv_opt(h, g_b), g_a, 'enable_disp', false);
f_rec_h_b_a = ft_rec_conv_opt(ft_rec_conv_opt(h, g_b), g_a, 'enable_disp', false);

%Plotting
figure;
tiledlayout(2,1);
nexttile
t = 1:numel(f_exp_h_a_b);
plot(t, f_exp_h_a_b, t, f_exp_h_b_a, t, f_exp_h_a_b - f_exp_h_b_a, 'LineWidth', 1.5);
grid on; axis tight;
legend({'$f(f(h, g_a), g_b)$', '$f(f(h, g_b), g_a)$', '$f(f(h, g_a), g_b) - f(f(h, g_b), g_a)$'}, ...
    'location', 'best', 'interpreter', 'latex');
title('Exp. Conv.');

nexttile
t = 1:numel(f_exp_h_a_b);
plot(t, f_rec_h_a_b, t, f_rec_h_b_a, t, f_rec_h_a_b - f_rec_h_b_a, 'LineWidth', 1.5);
grid on; axis tight;
legend({'$f(f(h, g_a), g_b)$', '$f(f(h, g_b), g_a)$', '$f(f(h, g_a), g_b) - f(f(h, g_b), g_a)$'}, ...
    'location', 'best', 'interpreter', 'latex');
title('Rec. Conv.');

%Error
err_exp_commute = err_SNMSE(f_exp_h_a_b, f_exp_h_b_a)
err_rec_commute = err_SNMSE(f_rec_h_a_b, f_rec_h_b_a)


