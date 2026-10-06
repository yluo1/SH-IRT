function plot_sh_fit_msq_omni_array_example
% Design multiple beamformers from directional frequency responses of pistons on spherical baffle oriented 
% over the horizontal plane to target a total omni-directional power response

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_odr = 12;
is_real = false;
oversample = 2;

% Generate directional frequency responses for spherical baffle piston oriented over azimuth on horizontal plane
N_dic = 4; % Number of dictionary elements
phi_dic = deg2rad(linspace(0, 360, N_dic + 1))'; % Equi-distributed steering angles over horizontal plane
phi_dic = phi_dic(1:end-1);
C_dic = zeros([(max_odr + 1)^2, N_dic]);
for n = 1:N_dic
    C_dic(:, n) = sh_enc_pist_sphere(max_odr, pi/2, phi_dic(n), 0.25, deg2rad(20), 400, 4);
end

% Constant power target
[theta, phi] = sh_grd_fib(oversample * (2 * max_odr + 1)^2); 
X_MS = ones(size(theta));

% Fit
[D_fit, C_fit, err, output] = sh_fit_msq(X_MS, theta, phi, 2 * max_odr, is_real, 'MOMS', 'B', C_dic);

beamformer_coeffs = output.W

dB_lim = [-20, 6];
h_dic = sh_plt(C_dic, 'mercator', is_real, 'dB_lim', [12, 20], 'title_name', 'Dictionary');
h_beam = sh_plt(C_fit, 'mercator', is_real, 'dB_lim', [-20, 6], 'title_name', 'Beam');
h_fit = sh_plt(D_fit, 'mercator', is_real, 'dB_lim', [-3 3], 'title_name_override', 'Total Power', 'disp_phase', false);
;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out_dir = 'figs/figs_sh/fit';
% if ~isfolder(out_dir)
%     mkdir(out_dir);
% end
% 
% for n = 1:N_dic
%     exportgraphics(h_dic{n}, fullfile(out_dir,  ['omni_array_dic_', num2str(n), '.png'] ));
%     exportgraphics(h_beam{n}, fullfile(out_dir,  ['omni_array_beam_', num2str(n), '.png'] ));
% end
% exportgraphics(h_fit{1}, fullfile(out_dir,  ['omni_array_fit.png'] ));
% 

