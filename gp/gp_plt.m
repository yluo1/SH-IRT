function h_fig_list = gp_plt(cov_name, options)
%Plot Gaussian process covariance functions

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%cov_name:          String, covariance function name {'cov_sqx_chw_ns', 'cov_sqx_chw_sqx'}

%options:           Struct

%options.cov_sqx_chw_ns_f_lo:       Scalar, frequency start  (Hz)
%options.cov_sqx_chw_ns_f_hi:       Scalar, Frequency end (Hz)
%options.cov_sqx_chw_ns_f_0:        [1 x N_f_0] Frequency f_0 (Hz)
%optons.cov_sqx_chw_ns_sigma:       Scalar, covariance scaling parameter
%optons.cov_sqx_chw_ns_ell_list:    [1 x N_ell], list of wavelength scaling
%optons.cov_sqx_chw_ns_gamma_list:  [1 x N_gamma], list of wavelength exponent

%options.enable_export:             Logical, if true, export figures to files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h_fig_list:        [1 x *] cell of figure handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Plot covariance function w.r.t. varying hyperparameters

%gp_plt('cov_sqx_chw_ns', 'enable_export', true);
%gp_plt('cov_sqx_chw_sqx', 'enable_export', true);

arguments
    cov_name (1,:) char {mustBeMember(cov_name, {'cov_sqx_chw_ns', 'cov_sqx_chw_sqx'})} = 'cov_sqx_chw_ns';

    options.cov_sqx_chw_ns_f_0_list  (1,:) double {mustBePositive} = [50, 500, 5000];
    options.cov_sqx_chw_ns_f_lo (1,1) double {mustBePositive} = 50;
    options.cov_sqx_chw_ns_f_hi (1,1) double {mustBePositive}  = 24000;
    options.cov_sqx_chw_ns_sigma (1,1) double {mustBeNonnegative} = 1;
    options.cov_sqx_chw_ns_ell_list (1,:) double {mustBePositive} = 343;
    options.cov_sqx_chw_ns_gamma_list (1,:) double {mustBePositive} =  [0.5, 1, 1.5];


    options.cov_sqx_chw_sqx_f_0_list  (1,:) double {mustBePositive} = [50, 500, 5000];
    options.cov_sqx_chw_sqx_f_lo (1,1) double {mustBePositive} = 50;
    options.cov_sqx_chw_sqx_f_hi (1,1) double {mustBePositive}  = 24000;
    options.cov_sqx_chw_sqx_sigma (1,1) double {mustBeNonnegative} = 1;
    options.cov_sqx_chw_sqx_ell_c_list (1,:) double {mustBePositive} = [1 2]/2;
    options.cov_sqx_chw_sqx_ell_f_list (1,:) double {mustBePositive} = [2 4]/2;


    options.fontsize (1,1) double {mustBePositive, mustBeInteger} = 20;
    options.FaceColor (1,:) char {mustBeMember(options.FaceColor, {'flat', 'interp'}) } = 'flat';
    options.dB_lim (1,2) double  = [-60 0];

    options.enable_export (1,1) logical = false;
end

if strcmp(cov_name, 'cov_sqx_chw_ns')

    N_freq = 600;
    N_theta_phi = 400; 

    f_0_list = options.cov_sqx_chw_ns_f_0_list;
    f_lo = options.cov_sqx_chw_ns_f_lo;
    f_hi = options.cov_sqx_chw_ns_f_hi;

    freq_list = logspace(log10(f_lo), log10(f_hi), N_freq);
    phi_list = linspace(0, pi, N_theta_phi);
    theta = pi / 2;

    [freq_mat, phi_mat] = meshgrid(freq_list, phi_list);
    X = [2 * pi * freq_mat(:), theta * ones(N_freq * N_theta_phi, 1), phi_mat(:)];


    %Hyperparameters
    sigma = options.cov_sqx_chw_ns_sigma;
    ell_list = options.cov_sqx_chw_ns_ell_list;
    gamma_list = options.cov_sqx_chw_ns_gamma_list;

    N_ell = numel(ell_list);
    N_gamma = numel(gamma_list);
    N_f_0 = numel(f_0_list);

    %Plotting
    h_fig = figure;
    h_fig.Position = [100, 100, 1600, 900];
    h_fig_list{1} = h_fig;

    %Setup tiled layouts
    if N_ell == 1 && N_f_0 == 1
        h_t = tiledlayout(2, ceil(N_gamma * N_f_0 /2));

    elseif N_ell == 1 && N_f_0 > 1
        h_t = tiledlayout(N_f_0, N_gamma);

    else
        h_t = tiledlayout(N_ell, N_gamma * N_f_0);
    end

    %Iterate over f0 
    for m = 1:N_f_0 

        X0 = [2 * pi * f_0_list(m), theta,  phi_list(1)];

        for i = 1:N_ell %Iterate over wavelength scale
            ell = ell_list(i);
    
            for j = 1:N_gamma %Iterate over wavelength exponent
                gamma = gamma_list(j);
                
                nexttile
                
                %Compute covariance
                K = cov_sqx_chw_ns(X, X0, sigma, ell, gamma);    
                K_phi_freq = reshape(K, [N_theta_phi, N_freq]);
    
                
                h_pc = pcolor(freq_list, rad2deg(phi_list), mag2db(K_phi_freq));
                set(h_pc, 'EdgeColor', 'none');
                set(h_pc, 'FaceColor', options.FaceColor);
                set(gca, 'XScale', 'log');
                %set(gca, 'YScale', 'log');
                %yticks([1, 10, 100]);

                h_cb = colorbar;
                ylabel(h_cb, 'dB', 'fontsize', options.fontsize);
                set(gca, 'fontsize', options.fontsize - 1);
            
                xlabel(['Frequency $f_1 \, | \, f_0$ = ', num2str(f_0_list(m)), ' (Hz)'], 'fontsize', options.fontsize, 'interpreter', 'latex');
                ylabel('Angle Dist. (Degrees)', 'fontsize', options.fontsize, 'interpreter', 'latex');
            
                title({['$\ell$ = ', num2str(ell), ', $\gamma$ = ', num2str(gamma)]}, ...
                    'fontsize', options.fontsize + 1, 'interpreter', 'latex');
    
                clim(options.dB_lim);
            end
        end
    end

    title(h_t, 'Squared Exponential Chordal Distance with Non-stationary Frequency', 'fontsize', options.fontsize + 1);


elseif strcmp(cov_name, 'cov_sqx_chw_sqx')

    N_freq = 600;
    N_theta_phi = 400; 

    f_0_list = options.cov_sqx_chw_sqx_f_0_list;
    f_lo = options.cov_sqx_chw_sqx_f_lo;
    f_hi = options.cov_sqx_chw_sqx_f_hi;

    freq_list = logspace(log10(f_lo), log10(f_hi), N_freq);

    phi_list = linspace(0, pi, N_theta_phi);
    theta = pi / 2;

    [freq_mat, phi_mat] = meshgrid(freq_list, phi_list);
    X = [2 * pi * freq_mat(:), theta * ones(N_freq * N_theta_phi, 1), phi_mat(:)];


    %Hyperparameters
    sigma = options.cov_sqx_chw_sqx_sigma;
    ell_c_list = options.cov_sqx_chw_sqx_ell_c_list;
    ell_f_list = options.cov_sqx_chw_sqx_ell_f_list;

    N_ell_c = numel(ell_c_list);
    N_ell_f = numel(ell_f_list);
    N_f_0 = numel(f_0_list);

    %Plotting
    h_fig = figure;
    h_fig.Position = [100, 100, 1600, 900];
    h_fig_list{1} = h_fig;

    %Setup tiled layouts
    if N_ell_c == 1 && N_f_0 == 1
        h_t = tiledlayout(2, ceil(N_ell_c * N_ell_f *  N_f_0 /2));
    else
        h_t = tiledlayout(N_f_0, N_ell_c * N_ell_f);
    end

    %Iterate over f0 
    for m = 1:N_f_0 

        X0 = [2 * pi * f_0_list(m), theta,  phi_list(1)];

        for i = 1:N_ell_c %Iterate over chordal distance scale
            ell_c = ell_c_list(i);
    
            for j = 1:N_ell_f %Iterate over log-frequency scale
                ell_f = ell_f_list(j);
                
                nexttile
                
                %Compute covariance
                K = cov_sqx_chw_sqx(X, X0, sigma, ell_c, ell_f);    
                K_phi_freq = reshape(K, [N_theta_phi, N_freq]);
    
                
                h_pc = pcolor(freq_list, rad2deg(phi_list), mag2db(K_phi_freq));
                set(h_pc, 'EdgeColor', 'none');
                set(h_pc, 'FaceColor', options.FaceColor);
                set(gca, 'XScale', 'log');
                %set(gca, 'YScale', 'log');
                %yticks([1, 10, 100]);

                h_cb = colorbar;
                ylabel(h_cb, 'dB', 'fontsize', options.fontsize);
                set(gca, 'fontsize', options.fontsize - 1);
            
                xlabel(['Frequency $f_1 \, | \, f_0$ = ', num2str(f_0_list(m)), ' (Hz)'], 'fontsize', options.fontsize, 'interpreter', 'latex');
                ylabel('Angle Dist. (Degrees)', 'fontsize', options.fontsize, 'interpreter', 'latex');
            
                title({['$\ell_c$ = ', num2str(ell_c), ', $\ell_f$ = ', num2str(ell_f)]}, ...
                    'fontsize', options.fontsize + 1, 'interpreter', 'latex');
    
                clim(options.dB_lim);
            end
        end
    end

    title(h_t, 'Squared Exponential Chordal Distance x Squared Exponential Log-Frequency Distance', 'fontsize', options.fontsize + 1);


else
    error('Unsupported cov_name');
end

%Export
if options.enable_export

    out_dir = 'figs/figs_t60';
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    for n = 1:numel(h_fig_list)
        exportgraphics(h_fig_list{n}, fullfile(out_dir, [cov_name, '_', num2str(n), '.png']));
    end
end