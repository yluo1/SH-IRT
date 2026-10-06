function plot_sh_fit_msq_null_array_example
% Design single beamformer from directional frequency responses of pistons on spherical baffle oriented 
% over the horizontal plane to target an omni-directional power response

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_odr = 12;
is_real = false;
oversample = 2;

% Target null pattern
theta_peak = pi/2;
phi_peak = deg2rad(-10);
ell = 0.5;

C_peak = sh_msq(sh_enc_rbf('SqExp', floor(max_odr / 2), theta_peak, phi_peak, ell, is_real), is_real);   
C_peak = C_peak / sh_dec(C_peak, theta_peak, phi_peak, is_real); %Normalize range be between [0, 1]
C_null = sh_enc_uni(max_odr) - C_peak; %1 - peak
D_tgt = sh_msq(C_null);

[theta, phi] = sh_grd_fib(oversample * (2 * max_odr + 1)^2); 
X_MS = sh_dec(D_tgt, theta, phi, is_real); % Power response target

% Generate directional frequency responses for spherical baffle piston oriented over azimuth on horizontal plane
N_dic = 12; % Number of dictionary elements
phi_dic = deg2rad(linspace(-180, 180, N_dic + 1))'; % Equi-distributed steering angles over horizontal plane
phi_dic = phi_dic(1:end-1);
C_dic = zeros([(max_odr + 1)^2, N_dic]);
for n = 1:N_dic
    C_dic(:, n) = sh_enc_pist_sphere(max_odr, pi/2, phi_dic(n), 0.25, deg2rad(20), 400, 4);
end

% Fit
K0 = (max_odr + 1)^2;
%idx_dic = [1, N_dic];
idx_dic = 1:N_dic;
N_idx_dic = numel(idx_dic);
[D_fit, C_fit, err, output] = sh_fit_msq(real(X_MS), theta, phi, 2 * max_odr, is_real, 'MP', ...
    'B', C_dic(:, idx_dic), 'B0', randn(N_idx_dic, K0) + randn(N_idx_dic, K0) * 1i);

h_tgt_pow = sh_plt(D_tgt, 'mercator', is_real, 'dB_lim', [-12, 6], 'disp_phase', false, 'title_name_override', 'Target Null Power');
h_fit_pow = sh_plt(D_fit, 'mercator', is_real, 'dB_lim', [-12, 6], 'disp_phase', false, 'title_name_override', 'Fitted Null Power');
h_fit = sh_plt(C_fit, 'mercator', is_real, 'dB_lim', [-12, 6], 'title_name', 'Fitted Null');

err = err_NSHMSQ(D_tgt, D_fit)

beamformer_coeffs = output.W

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out_dir = 'figs/figs_sh/fit';
% if ~isfolder(out_dir)
%     mkdir(out_dir);
% end
% 
% exportgraphics(h_tgt_pow{1}, fullfile(out_dir,  ['null_array_tgt_pow.png'] ));
% exportgraphics(h_fit_pow{1}, fullfile(out_dir,  ['null_array_fit_pow.png'] ));
% exportgraphics(h_fit{1}, fullfile(out_dir,  ['null_array_fit_resp.png'] ));
