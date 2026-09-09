function tst_com_matmul_exp_conv
%Test commutative properties of matrix multiplication and exponentiating convolution

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(623);
Fs = 48000;

%T = 1;
%T = 0.1;
T = 0.01;
M = ceil(T * Fs);

max_odr = 8;
N = (max_odr + 1)^2;

is_real = true;

B = randn((max_odr + 1)^2);
B2 = randn((max_odr + 1)^2);

C = sh_rand(max_odr, M, is_real);
g = ft_two_tap_FIR(0, -6, false);

%Fixed B
BC = B * C;
D_pre = zeros(N, M * numel(g)); %Pre-matrix-mult
for n = 1:N
    D_pre(n, :) = ft_exp_conv_opt(BC(n,:), g);
end

D_pst = zeros(N, M * numel(g)); %Post-matrix-mult
for n = 1:N
    D_pst(n, :) = ft_exp_conv_opt(C(n,:), g);
end
D_pst = B * D_pst;

err = err_SNMSE(D_pre, D_pst) %Does commute

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time-vary B
B2C = [B * C(:, 1:floor(M/2)), B2 * C(:, (floor(M/2)+1):end )];
D_pre = zeros(N, M * numel(g)); %Pre-matrix-mult
for n = 1:N
    D_pre(n, :) = ft_exp_conv_opt(B2C(n,:), g);
end

D_pst = zeros(N, M * numel(g)); %Post-matrix-mult
for n = 1:N
    D_pst(n, :) = ft_exp_conv_opt(C(n,:), g);
end
D_pst = [B * D_pst(:, 1:floor(M/2)), B2 * D_pst(:, (floor(M/2)+1):end )];

err = err_SNMSE(D_pre, D_pst) %Does not commute

figure; plot(D_pre(1, :)); hold on; plot(D_pst(1, :));