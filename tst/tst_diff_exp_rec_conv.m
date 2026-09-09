function tst_diff_exp_rec_conv
%Test differences in exponentiating and recursive convolutions

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(4123);
Fs = 48000;
T = 1;

M = ceil(T * Fs);
h = randn(1, M);

g = ft_two_tap_FIR(0, -6, false); %Larger differences

%g = ft_exp_design(6, 'RT60',  'enable_disp', true, 'Fs', 48000, 'RT60_sec', [2 0.35 0.4 0.2 0.3]);
%g = ft_exp_design(6, 'RT60',  'enable_disp', false, 'Fs', 48000, 'RT60_sec', [2 1e-4]);
%g = [0 g];


f_exp = ft_exp_conv_opt(h, g, 'enable_disp', true);
f_rec = ft_rec_conv_opt(h, g, 'enable_disp', true);

% f_exp = ft_exp_conv_direct(h, g);
% f_rec = ft_rec_conv_direct(h, g);

figure;
t = 1:numel(f_exp);
plot(t, f_exp, t, f_rec, t, f_exp - f_rec); grid on; axis tight;
legend({'exp', 'rec', 'exp - rec'}, 'location', 'best');

err = err_SNMSE(f_exp, f_rec) %Does not commute

