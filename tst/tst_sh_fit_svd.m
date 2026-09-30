function tst_sh_fit_svd(mode)
%Test sh_fit_svd.m conditions

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%mode:  String, test conditions {'picard_varycross', 'picard_domsingularvec'}
%           'picard_varycross':         Varying crossover tau
%           'picard_domsingularvec':    Varying concentration of first sigular vector
%                                       X_alpha = (1-alpha) * X + alpha * u_1 * sigma(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Test different conditions

%tst_sh_fit_svd('picard_varycross');
%tst_sh_fit_svd('picard_domsingularvec');

arguments
    mode (1,:) char {mustBeMember(mode, {'picard_varycross', 'picard_domsingularvec'})} = 'picard_varycross';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate and sample from random field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(441);
max_odr = 6;
is_real = false;
N_C = (max_odr + 1)^2;
N_pts = (max_odr + 1)^2;

%Reference field
C_ref = sh_rand(max_odr, 1, is_real);

%Sample at spherical coordinates
%[theta, phi] = sh_grd_fib(N_pts);       %Fibonnaci
[theta, phi] = sh_grd_rand(N_pts);  %Random over sphere

X = sh_dec(C_ref, theta, phi, is_real);
%Add noise
X = X + (randn(size(X)) + randn(size(X)) * 1i) * 5e-2;

%Plot reference
dB_lim = [-32, 24];

if strcmp(mode, 'picard_varycross')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Test mode = "picard" for varying crossover tau
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plot reference
    sh_plt(C_ref, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Reference', 'disp_theta_phi', [theta, phi]);
    [C_ls]  = sh_fit_svd(X, theta, phi, max_odr, is_real, 0, 'max', 'enable_disp', true, 'disp_dB_lim', dB_lim);
    
    %Picard cross criterion
    [C_pc_05]  = sh_fit_svd(X, theta, phi, max_odr, is_real, 0.5, 'picard', 'enable_disp', true, 'disp_dB_lim', dB_lim);
    [C_pc_025]  = sh_fit_svd(X, theta, phi, max_odr, is_real, 0.25, 'picard', 'enable_disp', true, 'disp_dB_lim', dB_lim);
    [C_pc_0125]  = sh_fit_svd(X, theta, phi, max_odr, is_real, 0.125, 'picard', 'enable_disp', true, 'disp_dB_lim', dB_lim);

elseif strcmp(mode, 'picard_domsingularvec')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Test robustness of mode = "picard" for varying concentrations of 
    % first sigular vector: X_alpha = (1-alpha) * X + alpha * u_1 * sigma(1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Y = sh_val(max_odr, theta, phi, is_real);  %[N x (max_odr + 1)^2]
    [U, S, V] = svd(Y, "econ");
    
    for alpha = [0, 0.5, 0.9, 0.99, 1]
        X_alpha = (1 - alpha) * X + alpha * U(:, 1) * S(1,1);
        sh_fit_svd(X_alpha, theta, phi, max_odr, is_real, 0.5, 'picard', 'enable_disp', true, 'disp_dB_lim', dB_lim);
    end

else
    error('Unsupported mode')
end
