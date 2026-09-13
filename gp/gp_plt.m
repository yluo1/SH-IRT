function h_fig_list = gp_plt(func_name, options)
%Plot Gaussian process prior functions

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%func_name:         String, prior mean or covariance function name
%                   {'mu_pow', 'mu_lpf', 'cov_sqx_chw_ns', 'cov_sqx_chw_sqx'}

%options:           Struct

% options.mu_pow_f_lo:                  Scalar, frequency start (Hz)
% options.mu_pow_f_hi:                  Scalar, Frequency end (Hz)
% options.mu_pow_alpha_list:            [1 x N_alpha], list of T60 at 0 Hz
% options.mu_pow_beta_list:             [1 x N_beta], list of decay exponents

% options.mu_lpf_f_lo:                  Scalar, frequency start (Hz)
% options.mu_lpf_f_hi:                  Scalar, Frequency end (Hz)
% options.mu_lpf_alpha_list:            [1 x N_alpha], list of T60 at 0 Hz
% options.mu_lpf_beta_list:             [1 x N_beta], list of decay exponents
% options.mu_lpf_omega_fc_list:         [1 x N_fc], list of cross-over frequencies

% options.cov_sqx_chw_ns_f_0_list:       [1 x N_f_0] Frequency f_0 (Hz)
% options.cov_sqx_chw_ns_f_lo:           Scalar, frequency start (Hz)
% options.cov_sqx_chw_ns_f_hi:           Scalar, Frequency end (Hz)
% options.cov_sqx_chw_ns_sigma:          Scalar, covariance scaling parameter
% options.cov_sqx_chw_ns_ell_list:       [1 x N_ell], list of wavelength scaling
% options.cov_sqx_chw_ns_gamma_list:     [1 x N_gamma], list of wavelength exponents
% 
% options.cov_sqx_chw_sqx_f_0_list:      [1 x N_f_0] Frequency f_0 (Hz)
% options.cov_sqx_chw_sqx_f_lo:          Scalar, frequency start (Hz)
% options.cov_sqx_chw_sqx_f_hi:          Scalar, Frequency end (Hz)
% options.cov_sqx_chw_sqx_sigma:         Scalar, covariance scaling parameter
% options.cov_sqx_chw_sqx_ell_c_list:    [1 x N_ell], list of wavelength scaling
% options.cov_sqx_chw_sqx_ell_f_list:    [1 x N_gamma], list of wavelength exponents
% 
% options.enable_export:                 Logical, if true, export figures to files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h_fig_list:        [1 x *] cell of figure handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Plot function w.r.t. varying hyperparameters

%gp_plt('mu_pow', 'enable_export', true);
%gp_plt('mu_lpf', 'enable_export', true);

%gp_plt('cov_sqx_chw_ns', 'enable_export', true);
%gp_plt('cov_sqx_chw_sqx', 'enable_export', true);

arguments
    func_name (1,:) char {mustBeMember(func_name, {'mu_pow', 'mu_lpf', 'cov_sqx_chw_ns', 'cov_sqx_chw_sqx'})} = 'cov_sqx_chw_ns';

    options.mu_pow_f_lo (1,1) double {mustBePositive} = 50;
    options.mu_pow_f_hi (1,1) double {mustBePositive}  = 24000;
    options.mu_pow_alpha_list (1,:) double {mustBeNonnegative} = [0.25, 1, 1.5];
    options.mu_pow_beta_list (1,:) double = [0, 0.125, 0.25];

    options.mu_lpf_f_lo (1,1) double {mustBePositive} = 50;
    options.mu_lpf_f_hi (1,1) double {mustBePositive}  = 24000;
    options.mu_lpf_alpha_list (1,:) double {mustBeNonnegative} = [0.25, 1, 1.5];
    options.mu_lpf_beta_list (1,:) double  {mustBeNonnegative} = [0, 0.25, 0.5];
    options.mu_lpf_omega_fc_list (1,:) double = 2 * pi * [50, 500, 5000];


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

if strcmp(func_name, 'mu_pow')

    N_freq = 32;
    f_lo = options.mu_pow_f_lo;
    f_hi = options.mu_pow_f_hi;
    freq_list = logspace(log10(f_lo), log10(f_hi), N_freq);
    omega = 2 * pi * freq_list;
       
    fontsize = options.fontsize - 2;

    alpha_list = options.mu_pow_alpha_list;
    beta_list = options.mu_pow_beta_list;

    N_alpha = numel(alpha_list);
    N_beta = numel(beta_list);

    mu_list = zeros([N_freq, N_alpha, N_beta]);
    for i = 1:N_alpha
        for j = 1:N_beta
            mu_list(:, i, j) = mu_pow(omega, alpha_list(i), beta_list(j));
        end
    end
    
    %Plot
    h_fig = figure;
    h_fig.Position = [100, 100, 900, 450];
    h_fig_list{1} = h_fig;

    line_style_list = {'-', '--', '-.', ':'};
    marker_style_list = {'*', 'o', 's', 'd'};
    color_list = orderedcolors("gem");
    legend_str = {};

    for i = 1:N_alpha
        for j = 1:N_beta
            semilogx(freq_list, mu_list(:, i, j), ...
                'Color', color_list(mod(i-1, size(color_list, 1)) + 1, :), ...
                'LineStyle', line_style_list{mod(j-1, numel(line_style_list) ) + 1 }, ...
                'Marker', marker_style_list{mod(j-1, numel(marker_style_list) ) + 1 }, 'MarkerSize', 8, ...
                'linewidth', 1.5);
            hold on;
            legend_str{end+1} = ['$\alpha = $ ', num2str(alpha_list(i)), ', $\beta = $ ', num2str(beta_list(j))];
        end
    end
    ylim([0, max(mu_list(:))]);
    grid on;
    axis tight;
    xlabel('Frequency (Hz)', 'fontsize', fontsize);
    ylabel('T60 (Seconds)', 'fontsize', fontsize);
    title('Prior Mean Power Function', 'fontsize', fontsize + 1);
    set(gca, 'fontsize', fontsize - 1);
    h_lg = legend(legend_str, 'location', 'best', 'NumColumns', N_alpha, 'interpreter', 'latex');
    set(h_lg, 'fontsize', fontsize - 1);
 
elseif strcmp(func_name, 'mu_lpf')

    N_freq = 32;
    f_lo = options.mu_lpf_f_lo;
    f_hi = options.mu_lpf_f_hi;
    freq_list = logspace(log10(f_lo), log10(f_hi), N_freq);
    omega = 2 * pi * freq_list;
       
    fontsize = options.fontsize - 2;

    alpha_list = options.mu_lpf_alpha_list;
    beta_list = options.mu_lpf_beta_list;
    omega_fc_list =  options.mu_lpf_omega_fc_list;

    N_alpha = numel(alpha_list);
    N_beta = numel(beta_list);
    N_fc = numel(omega_fc_list);

    mu_list = zeros([N_freq, N_alpha, N_beta, N_fc]);
    for i = 1:N_alpha
        for j = 1:N_beta
            for k = 1:N_fc
                mu_list(:, i, j, k) = mu_lpf(omega, alpha_list(i), beta_list(j), omega_fc_list(k));
            end
        end
    end
    
    %Plot
    h_fig = figure;
    h_fig.Position = [100, 100, 1600, 900];
    h_fig_list{1} = h_fig;

    line_style_list = {'-', '--', '-.', ':'};
    marker_style_list = {'*', 'o', 's', 'd'};
    color_list = orderedcolors("gem");
    legend_str = {};

    h_t = tiledlayout(1, N_fc);
    for k = 1:1:N_fc
        nexttile;
        for i = 1:N_alpha
            for j = 1:N_beta
                semilogx(freq_list, mu_list(:, i, j, k), ...
                    'Color', color_list(mod(i-1, size(color_list, 1)) + 1, :), ...
                    'LineStyle', line_style_list{mod(j-1, numel(line_style_list) ) + 1 }, ...
                    'Marker', marker_style_list{mod(j-1, numel(marker_style_list) ) + 1 }, 'MarkerSize', 8, ...
                    'linewidth', 1.5);
                hold on;
                legend_str{end+1} = ['$\alpha = $ ', num2str(alpha_list(i)), ', $\beta = $ ', num2str(beta_list(j))];
            end
        end

        ylim([0, max(mu_list(:))]);
        grid on;
        axis tight;
        xlabel('Frequency (Hz)', 'fontsize', fontsize);
        ylabel('T60 (Seconds)', 'fontsize', fontsize);
        title(['$f_c = $ ', num2str(omega_fc_list(k) / (2 * pi) ), ' Hz'], 'fontsize', fontsize + 1, 'interpreter', 'latex');
        set(gca, 'fontsize', fontsize - 1);
        h_lg = legend(legend_str, 'location', 'best', 'NumColumns', 1, 'interpreter', 'latex');
        set(h_lg, 'fontsize', fontsize - 1);
     
    end

    title(h_t, 'Prior Mean Low-Pass Function', 'fontsize', fontsize + 1);

elseif strcmp(func_name, 'cov_sqx_chw_ns')

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

    title(h_t, 'Squared Exponential Chordal Distance with Non-stationary Frequency Covariance Function', 'fontsize', options.fontsize + 1);


elseif strcmp(func_name, 'cov_sqx_chw_sqx')

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

    title(h_t, 'Squared Exponential Chordal Distance x Squared Exponential Log-Frequency Distance Covariance Function', 'fontsize', options.fontsize + 1);


else
    error('Unsupported func_name');
end

%Export
if options.enable_export

    out_dir = 'figs/figs_t60';
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    for n = 1:numel(h_fig_list)
        exportgraphics(h_fig_list{n}, fullfile(out_dir, [func_name, '_', num2str(n), '.png']));
    end
end