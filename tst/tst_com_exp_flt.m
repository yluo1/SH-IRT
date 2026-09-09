function tst_com_exp_flt
%Test commutative properties of convolution and exponentiating convolution

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(4123);
Fs = 48000;

%T = 1;
%T = 0.1;
T = 0.01;

M = ceil(T * Fs);
h = randn(1, M);

g = ft_two_tap_FIR(-1e-3, -3, false);   %Larger differences

u = randn(1, 10);

f_pre = ft_exp_conv_opt(conv(u, h), g);
f_pst = conv(u, ft_exp_conv_opt(h, g));
[f_pre, f_pst] = trunc_mat_min(f_pre, f_pst);

figure;
t = 1:numel(f_pre);
plot(t, f_pre, t, f_pst, t, f_pre - f_pst, 'linewidth', 1.5); 
grid on; axis tight;
legend({'pre', 'pst', 'pre - pst'}, 'location', 'best');

err = err_SNMSE(f_pre, f_pst)
