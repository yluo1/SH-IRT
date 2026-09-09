function tst_sh_pdf_est
%Test density fititng methods for sh_pdf_est.m

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(11213); %original

%rng(113);
%rng(117);
%rng(1178);

%max_odr = 8;
max_odr = 6; %Original
%max_odr = 4;
%max_odr = 3;
%max_odr = 2;
%max_odr = 1;

is_real = true;

%N = 10;
%N = 50;
N = 100; %Default
%N = 250;

theta = acos(2 * rand(N, 1) - 1); 
phi = rand(N, 1) * 2 * pi;
%dB_lim = [-60, 0];
dB_lim = [-48, 0];
fontsize = 20;

%theta = [pi/2; pi/2; pi/2; pi/2];
%phi = [0; pi/4; pi/8;  pi];

% theta = acos(2 * rand(N, 1) - 1) / 2; %[0 to pi/2]
% phi = rand(N, 1) * 2 * pi;
% theta = [theta; pi]; %Add south pole point
% phi = [phi; 0];


% mode_list = {'Uniform', 'KernelDensityProjPost', 'KernelDensitySqProj', 'PrincipalAxes', 'MaxLogLike'};
% mode_list_name = {'Uniform', 'Kernel Density Projection', 'Kernel Density Squared Projection', 'Principal Axes', 'Maximum Log-Likelihood'};

% mode_list = {'Uniform', 'KernelDensityProjPost', 'KernelDensitySqProj', 'MaxLogLike'};
% mode_list_name = {'Uniform', 'Projection Kernel', 'Squared Projection Kernel', 'Maximum Log-Likelihood'};

mode_list = {'Uniform', 'KernelDensityProjPost', 'SpectralMaxLogLike', 'MaxLogLike'};
mode_list_name = {'Uniform', 'Projection Kernel', 'Maximum Log-Spectral', 'Maximum Log-Likelihood'};


N_modes = numel(mode_list);
C_list = cell(1, N_modes);
D_list = cell(1, N_modes);
loglike_list = cell(1, N_modes);
fig_list = cell(1, N_modes);

N_x0 = 100;
if is_real
    MaxLogLike_x0_list = randn((max_odr + 1)^2, N_x0);
else
    MaxLogLike_x0_list = randn((max_odr + 1)^2, N_x0) + randn((max_odr + 1)^2, N_x0) * 1i;
end
%MaxLogLike_x0_list = [];

MaxLogLike_relax_eq_constr = false;
%MaxLogLike_relax_eq_constr = true;

%MaxLogLike_mode = 'fmincon';
MaxLogLike_mode = 'sdp';

for i = 1:N_modes

    if strcmp(mode_list{i}, 'MaxLogLike')
%        [C_list{i}, D_list{i}, loglike_list{i}] = sh_pdf_est(theta, phi, max_odr, is_real, mode_list{i}, 'MaxLogLike_x0_list', MaxLogLike_x0_list);
        
        %Include all previous solutions as starting point
        [C_list{i}, D_list{i}, loglike_list{i}] = sh_pdf_est(theta, phi, max_odr, is_real, mode_list{i}, ...
            'MaxLogLike_mode', MaxLogLike_mode, ...
            'MaxLogLike_x0_list', [cell2mat(C_list(1:i-1)), MaxLogLike_x0_list], ...
            'MaxLogLike_relax_eq_constr', MaxLogLike_relax_eq_constr);

    else
        [C_list{i}, D_list{i}, loglike_list{i}] = sh_pdf_est(theta, phi, max_odr, is_real, mode_list{i});
    end
    fig_list{i} = sh_plt(D_list{i}, 'mercator', is_real, ...
        'disp_phase', false, 'disp_theta_phi', [theta, phi], 'dB_lim', dB_lim, 'title_name', mode_list_name{i}, 'fontsize', fontsize, ...
        'fig_size', [560, 320]);
    loglike_list{i}
    title([mode_list_name{i}, ' Density'], 'fontsize', fontsize);
end
;

%Compute CDF
N_fwd = 100;
theta_fwd = linspace(0, pi, N_fwd);
N_inv = 11;
CDF_C = sh_cdf_theta(D_list{end}, theta_fwd);
u = linspace(0, 1, N_inv)';
theta_CDF_inv = sh_cdf_inv_theta(D_list{end}, u);

%Plotting
fontsize = 26;
h_cdf = figure;
h_cdf.Position = [100, 100, [500, 600] * 4/5];
plot(rad2deg(theta_fwd), CDF_C, '-', rad2deg(theta_CDF_inv), u, '*', 'linewidth', 4, 'markersize', 14); grid on; axis tight; 
h_lg = legend('$u = F_{\Theta}(\theta)$', '$\theta_* = F_{\Theta}^{-1}(u_*)$', ...
    'location', 'best', 'interpreter', 'latex');
set(h_lg, 'fontsize', fontsize );
xlabel('Co-latitude $\theta$ (Degree)', 'fontsize', fontsize, 'interpreter', 'latex');
ylabel('CDF $F_{\Theta}(\theta)$', 'fontsize', fontsize, 'interpreter', 'latex');
set(gca, 'fontsize', fontsize - 1);
title('Cumulative Distribution', 'fontsize', fontsize + 1, 'interpreter', 'latex');


out_dir = 'figs_sh';
if ~isfolder(out_dir)
    mkdir(out_dir);
end

exportgraphics(h_cdf, fullfile(out_dir, ['cdf_', mode_list{end}, '.png']));

for i = 1:N_modes
    exportgraphics(fig_list{i}{1}, fullfile(out_dir, [mode_list{i}, '.png']));
end