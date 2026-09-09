function plot_sh_mul_example
%Compare function expansion with SH multiplication operator's singular value exponentiation

P_C = 3;
P_D = 3;

% P_C = 5;
% P_D = 5;

theta = pi/2;
phi = 0;

%ell = 1;
ell = 0.5;

C = sh_enc_rbf('SqExp', P_C, theta, phi, ell, false, false);
D = sh_enc_uni(P_D);

[E, A] = sh_mul(C, D); 

B = A(1:size(A, 2), :); %Truncated A
B = (B + B')/2;         %Ensure Hermitian
[eig_V, eig_D] = svd(B);
err_eig = norm(eig_V * eig_D * eig_V' - B)

[U, S, V] = svd(A);
err_svd = norm(U * S * V' - A)

%sh_plt(U * S * V' * D)

%Exponentiation equivalence
%alpha = 1;     %Identity

%alpha = 0.75;    %Broaden
alpha = 0.5;    %Broaden
%alpha = 0.005;    %Broaden

%alpha = 2;         %Tighten
%alpha = 2.5;       %Tighten
%alpha = 4;         %Tighten
%alpha = 8;         %Tighten
%alpha = 16;         %Tighten

C_alpha = sh_enc_rbf('SqExp', P_C, theta, phi, ell / sqrt(alpha), false, false);

%E_alpha = U * S.^(alpha) * V' * D; %From singular value decomposition
A_alpha = eig_V * eig_D.^(alpha) * eig_V' * D; %From Eigendecomposition


% dB_lim = [-6, 6];
% sh_plt(C_alpha, 'mercator', false, 'title_name', 'C alpha', 'dB_lim', dB_lim)
% sh_plt(E_alpha, 'mercator', false, 'title_name', 'Exp', 'dB_lim', dB_lim);

N_val = 256;
theta_val = pi/2 * ones(N_val, 1);
phi_val = linspace(0, pi, N_val)';
f_C_alpha = real(sh_dec(C_alpha, theta_val, phi_val));
f_A_alpha = real(sh_dec(A_alpha, theta_val, phi_val));
f_ref_alpha = rbf_val('SqExp', 2 * sin(abs(phi_val)/2), ell / sqrt(alpha), false);

err_C_A_alpha = norm(f_C_alpha - f_A_alpha)
err_C_ref_alpha = norm(f_C_alpha - f_ref_alpha)
err_A_ref_alpha = norm(f_A_alpha - f_ref_alpha)


fontsize = 14;
figure;
plot(phi_val, f_C_alpha, phi_val, f_A_alpha, phi_val, f_ref_alpha, 'k',  'linewidth', 1.5);
xlabel('d', 'fontsize', fontsize);
ylabel('f(d)', 'fontsize', fontsize);
grid on; axis tight;
h_l = legend('Re-expansion', 'Exponentiation', 'Reference');
set(h_l, 'fontsize', fontsize - 1);
set(gca, 'fontsize', fontsize - 1);

% %Exponentiate the singular values
% plot_sph_coeffs((eig_V * eig_D * eig_V') * D, [], 'all', [], [], 'disp_mode', 'mag', 'dB_lim', [-30, 10]);
% plot_sph_coeffs((eig_V * (eig_D.^(0.5)) * eig_V') * D, [], 'all', [], [], 'disp_mode', 'mag', 'dB_lim', [-30, 10]);
% plot_sph_coeffs((eig_V * (eig_D.^(0.25)) * eig_V') * D, [], 'all', [], [], 'disp_mode', 'mag', 'dB_lim', [-30, 10]);
% plot_sph_coeffs((eig_V * (eig_D.^(2)) * eig_V') * D, [], 'all', [], [], 'disp_mode', 'mag', 'dB_lim', [-30, 10]);
% plot_sph_coeffs((eig_V * (eig_D.^(4)) * eig_V') * D, [], 'all', [], [], 'disp_mode', 'mag', 'dB_lim', [-30, 10]);