function tst_bf_mvdr
% Test bf_mvdr.m condition number modes

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(41735);
N = 10;
v = randn(N, 1) + randn(N, 1)*1i;
R = randn(N) + randn(N)*1i; R = R*R';
R = R - 0.999 * eigs(R, 1, 'smallestreal') * eye(N);
%[v_max, d_max] = eigs(R, 1, 'largestreal');     R = R + 999 * v_max * (v_max');
N_kappa = 64;
kappa_list = logspace(log10(1.1), log10(100000), N_kappa);
W_diag = zeros(N, N_kappa);   SINR_diag = zeros(1, N_kappa);     G_diag = zeros(1, N_kappa);
W_thres = zeros(N, N_kappa);  SINR_thres = zeros(1, N_kappa);    G_thres = zeros(1, N_kappa);
W_exp = zeros(N, N_kappa);    SINR_exp = zeros(1, N_kappa);      G_exp = zeros(1, N_kappa);
W_sse = zeros(N, N_kappa);    SINR_sse = zeros(1, N_kappa);      G_sse = zeros(1, N_kappa);
W_ase = zeros(N, N_kappa);    SINR_ase = zeros(1, N_kappa);      G_ase = zeros(1, N_kappa);

for n = 1:N_kappa
    [W_diag(:, n), ~, SINR_diag(n), G_diag(n)]      = bf_mvdr(v, R, kappa_list(n), 'diag');
    [W_thres(:, n), ~, SINR_thres(n), G_thres(n)]   = bf_mvdr(v, R, kappa_list(n), 'thres');
    [W_exp(:, n), ~, SINR_exp(n), G_exp(n)]         = bf_mvdr(v, R, kappa_list(n), 'exp');
    [W_sse(:, n), ~, SINR_sse(n), G_sse(n)]         = bf_mvdr(v, R, kappa_list(n), 'sse');
    [W_ase(:, n), ~, SINR_ase(n), G_ase(n)]         = bf_mvdr(v, R, kappa_list(n), 'ase');
end

%Plotting
fontsize = 14;
h_fig = figure;
h_fig.Position = [100, 100, 800, 600];
tiledlayout(2,1);
nexttile;
semilogx(kappa_list, pow2db(G_diag), kappa_list, pow2db(G_thres), kappa_list, pow2db(G_exp), kappa_list, pow2db(G_sse),  kappa_list, pow2db(G_ase),  'linewidth', 1.5); 
grid on; axis tight;
xlabel('Max Condition Number \kappa Constraint', 'fontsize', fontsize);
ylabel('Electrical Gain (dB)', 'fontsize', fontsize);
title('MVDR Electrical Gain across Conditioning Methods', 'fontsize', fontsize + 1); set(gca, 'fontsize', fontsize - 1);
h_lg = legend('diag', 'thres', 'exp', 'sse', 'ase', 'location', 'best'); set(h_lg, 'fontsize', fontsize - 1);

nexttile;
semilogx(kappa_list, pow2db(SINR_diag), kappa_list, pow2db(SINR_thres), kappa_list, pow2db(SINR_exp), kappa_list, pow2db(SINR_sse), kappa_list, pow2db(SINR_ase), 'linewidth', 1.5); 
grid on; axis tight;
xlabel('Max Condition Number \kappa Constraint', 'fontsize', fontsize);
ylabel('Directivity Index (dB)', 'fontsize', fontsize);
title('MVDR Signal-to-interference-plus-noise across Conditioning Methods', 'fontsize', fontsize + 1); set(gca, 'fontsize', fontsize - 1);
h_lg = legend('diag', 'thres', 'exp', 'sse', 'ase', 'location', 'best'); set(h_lg, 'fontsize', fontsize - 1);