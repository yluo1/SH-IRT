function [log_T60, h_fig] = gp_t60_sample(omega, theta, phi, num_evals, obs, options)
%Sample log T60(omega, theta, phi) functions over an evaluation grid from a Gaussian process

%Evaluation grid is Cartesian product of {omega x {theta, phi}} angular frequency x spherical coordinates

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%omega:                [N_B x 1]    Angular frequency
%theta:                [N_E x 1]    Co-latitude [0, pi]
%phi:                  [N_E x 1]    Azimuth [0, 2 * pi)

%num_evals:             Number of times to sample evaluation grid

%obs:                   Struct, observed T60(omega, theta, phi), see gp_obs_opts.m
%                       [] for none

%options:               Struct, 
%options.options_mu:    Struct, prior mean options, see gp_mu_opts.m
%options.options_cov:   Struct, prior covariance options, see gp_cov_opts.m

%options.sample_method: String, sampling method {'mvnrnd', 'mean'}
%                       'mvnrnd':   Sample from multivariate gausian
%                       'mean'      Sample the mean

%options.enable_disp:   Logical, if true, display plots
%options.options_disp:  Struct, display options, see gp_disp_opts.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%log_T60:           [N_B x N_E x num_evals]
%h_fig:             Handle to figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Draw IID sample log-T60 functions from GP prior at common spherical coordinate
    
% options_mu    = gp_mu_opts('mu_func', 'power', 'mu_alpha', 1, 'mu_beta', 0.25, 'mu_fc', 8000);
% options_cov   = gp_cov_opts('cov_sigma', sqrt(2)/2, 'cov_gamma', 2/3, 'cov_ell', 343 * 1);

% N_B = 1024;
% N_E = 1;

% freq          = logspace(log10(20), log10(24000), N_B)';
% omega         = 2 * pi * freq;
% [theta, phi]  = sh_fib(N_E);

% rng(21136);
% num_evals = 5;
% [log_T60_prior] = gp_t60_sample(omega, theta, phi, num_evals, [], 'enable_disp', true, ...
%   'options_mu', options_mu, 'options_cov', options_cov);

arguments
    omega (:,1) double {mustBeNonnegative} = linspace(0, pi, 16)';
    theta (:,1) double = [0];
    phi   (:,1) double = [0];

    num_evals (1,1) double {mustBeNonnegative} = 1;

    %Observations
    obs = [];

    %Priors
    options.options_mu  = gp_mu_opts;
    options.options_cov = gp_cov_opts;

    options.sample_method (1,:) char {mustBeMember(options.sample_method, {'mvnrnd', 'mean'})} = 'mvnrnd';

    %Display
    options.enable_disp (1,1) logical = false;
    options.options_disp = gp_disp_opts;

end

N_B = numel(omega);
N_E = numel(theta);
N_BE = N_B * N_E;

assert(numel(theta) == numel(phi), 'theta, phi size mismatch');

%Compute evaluation grid
freq = max(1, omega / (2 * pi)); %[N_B x 1]
freq_grid = repmat(freq, [1, N_E]);
theta_grid  = repmat(theta', [N_B, 1]);
phi_grid    = repmat(phi', [N_B, 1]);
X = [2 * pi * freq_grid(:), theta_grid(:), phi_grid(:)]; %[N_B * N_E x 1]

%Check observations
has_obs = ~isempty(obs);
if has_obs 
    has_obs = numel(obs.theta) > 0;
end

%Compute prior mean on evaluation grid
mu_omega = gp_mu(omega, options.options_mu);
mu_prior = log(mu_omega);
mu = repmat(mu_prior, [N_E, 1]); %[N_B * N_E x 1]

%Compute prior covariance on evaluation grid
if strcmp(options.sample_method, 'mvnrnd') %Sample from multivariate gaussian

    Sigma = gp_cov(X, X, options.options_cov); %[N_B * N_E x N_B * N_E] full covariance

    %Compute observation covariance matrix inverse
    if has_obs
        N_S = numel(obs.omega);
        assert(all([numel(obs.omega), numel(obs.theta), numel(obs.phi), numel(obs.T60), numel(obs.log_noise_std)] == N_S), 'obs fields mismatch');
    
        freq_obs = max(1, obs.omega / (2 * pi) ); %[N_S x 1]
        X_obs = [2 * pi * freq_obs(:), obs.theta(:), obs.phi(:)];
        
        Sigma_obs = gp_cov(X_obs, X_obs, options.options_cov); %[N_S x N_S]
        KXobs     = gp_cov(X, X_obs, options.options_cov);  %[N_B * N_E x N_S]
    
        %Add noise variance
        Sigma_obs = Sigma_obs + diag(obs.log_noise_std.^2);
    
        %Truncated 0 eigenvalues, and invert
        [eig_V_obs, eig_D_obs] = eig(Sigma_obs);
        eig_list_obs = diag(eig_D_obs);    
        % min_eig_val_obs = min(eig_list_obs)
        % max_eig_val_obs = max(eig_list_obs)
        eig_list_obs(eig_list_obs > 0) = 1 ./ eig_list_obs(eig_list_obs > 0);
        eig_list_obs(eig_list_obs <= 0) = 0;
        eig_D_obs = diag(eig_list_obs);    
        Sigma_inv_obs = eig_V_obs * eig_D_obs * eig_V_obs'; %[N_S x N_S]
    
        %Compute mean prior
        mu_omega_obs = gp_mu(obs.omega, options.options_mu);
        y = log(obs.T60) - log(mu_omega_obs);
        
        %Update mean
        KXobsSigma = KXobs * Sigma_inv_obs; %[N_B * N_E x N_S] intermediary
        mu = mu + KXobsSigma * y;
    
        %Update covariance        
        Sigma = Sigma - KXobsSigma * KXobs';
    
    end
    
    %Truncated 0 eigenvalues
    [eig_V, eig_D] = eig(Sigma);
    eig_list = diag(eig_D);
    % min_eig_val = min(eig_list)
    % max_eig_val = max(eig_list)
    eig_list(eig_list < 0) = 0;
    eig_D = diag(eig_list);
    
    Sigma = real(eig_V * eig_D * eig_V');
    Sigma = Sigma + options.options_cov.cov_jit * eye(size(Sigma));
    Sigma = (Sigma + Sigma')/2;     %Force symmetric
    
    %Sample GP
    log_T60 = reshape( mvnrnd(mu, Sigma, num_evals)', [N_B, N_E, num_evals]); % [N_B x N_E x num_evals]    
    var_log_T60 = diag(Sigma);


elseif strcmp(options.sample_method, 'mean')   %Sample from mean

    Sigma_diag = zeros(N_BE, 1);               %[N_B * N_E x 1] diagonal covariance
    for n = 1:N_BE
        Sigma_diag(n) = gp_cov(X(n), X(n), options.options_cov);
    end

    %Compute observation covariance matrix inverse
    if has_obs
        N_S = numel(obs.omega);
        assert(all([numel(obs.omega), numel(obs.theta), numel(obs.phi), numel(obs.T60), numel(obs.log_noise_std)] == N_S), 'obs fields mismatch');
    
        freq_obs = max(1, obs.omega / (2 * pi) ); %[N_S x 1]
        X_obs = [2 * pi * freq_obs(:), obs.theta(:), obs.phi(:)];
        
        Sigma_obs = gp_cov(X_obs, X_obs, options.options_cov); %[N_S x N_S]
        KXobs     = gp_cov(X, X_obs, options.options_cov);  %[N_B * N_E x N_S]
    
        %Add noise variance
        Sigma_obs = Sigma_obs + diag(obs.log_noise_std.^2);
    
        %Truncated 0 eigenvalues, and invert
        [eig_V_obs, eig_D_obs] = eig(Sigma_obs);
        eig_list_obs = diag(eig_D_obs);    
        % min_eig_val_obs = min(eig_list_obs)
        % max_eig_val_obs = max(eig_list_obs)
        eig_list_obs(eig_list_obs > 0) = 1 ./ eig_list_obs(eig_list_obs > 0);
        eig_list_obs(eig_list_obs <= 0) = 0;
        eig_D_obs = diag(eig_list_obs);    
        Sigma_inv_obs = eig_V_obs * eig_D_obs * eig_V_obs'; %[N_S x N_S]
    
        %Compute mean prior
        mu_omega_obs = gp_mu(obs.omega, options.options_mu);
        y = log(obs.T60) - log(mu_omega_obs);
        
        %Update mean
        KXobsSigma = KXobs * Sigma_inv_obs; %[N_B * N_E x N_S] intermediary
        mu = mu + KXobsSigma * y;
    
        %Update covariance
        Sigma_diag = Sigma_diag - sum(KXobsSigma .* KXobs, 2);
    
    end
    
    %Truncated 0 eigenvalues
    % min_eig_val = min(Sigma_diag)
    % max_eig_val = max(Sigma_diag)
    Sigma_diag(Sigma_diag < 0) = 0;
    Sigma_diag = Sigma_diag + options.options_cov.cov_jit;

    %Sample GP mean
    log_T60 = repmat(reshape(mu, [N_B, N_E]), [1, 1, num_evals]);
    var_log_T60 = Sigma_diag;

else
    error('Unsupported options.sample_method');
end
%Undo log-transform
mu_T60 = exp(mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_fig = [];
if options.enable_disp    
    if ~has_obs
        name_mod = 'Prior';
        disp_mu_T60 = mu_T60(1:N_B);
        disp_var_log_T60 = var_log_T60(1:N_B);
    else
        name_mod = 'Posterior';
        disp_mu_T60 = reshape(mu_T60, [N_B, N_E]);
        disp_var_log_T60 = reshape(var_log_T60, [N_B, N_E]);
    end

    fontsize = options.options_disp.disp_fontsize;

    h_fig = figure;
    h_fig.Position = options.options_disp.disp_position;
    legend_str = {};

    colororder(options.options_disp.disp_colororder); %Set color order
    xscale('log'); hold on;

    %Plot prior/posterior mean and variance
    if N_E == 1 || ~has_obs %Only one evaluation angle, or is prior
        
        if options.options_disp.disp_mean
            semilogx(freq, disp_mu_T60, 'k--', 'linewidth', 3); hold on; 
            legend_str{end+1} = ['Prior Mean'];
        end

        if options.options_disp.disp_var
            fill([freq; flipud(freq)]', [exp( (log(disp_mu_T60) + 2 * sqrt(disp_var_log_T60) ));  flipud(exp( (log(disp_mu_T60) - 2 * sqrt(disp_var_log_T60) ))) ]',  ...
                [0.5 0.5 0.5], 'FaceAlpha', options.options_disp.disp_var_transparency, 'EdgeAlpha', 0); hold on;
            legend_str{end+1} = [name_mod, ' Mean $\pm$ $2$ Std.'];
        end

    else %multiple evaluation angles and is posterior

        for m = 1:options.options_disp.disp_sample_stride:N_E

            if options.options_disp.disp_mean
                if options.options_disp.disp_legend_mean
                    vis_mn = 'on';
                    if options.options_disp.disp_legend_compact
                        legend_str{end+1} = ['$(', num2str(rad2deg(theta(m)), 3), '^\circ,' num2str(rad2deg(phi(m)), 3), '^\circ)$'];    
                    else
                        legend_str{end+1} = ['Posterior Mean ', '$T_{60}(\omega,', num2str(rad2deg(theta(m)), 3), '^\circ,' num2str(rad2deg(phi(m)), 3), '^\circ)$'];
                    end
                else %options.options_disp.disp_legend_mean is false, display one entry only
                    if m == 1
                        legend_str{end+1} = ['Posterior Mean'];
                        vis_mn = 'on';
                    else
                        vis_mn = 'off';
                    end
                end
                semilogx(freq, disp_mu_T60(:, m), '--', 'linewidth', 3, 'HandleVisibility', vis_mn); hold on; 
            end

            if options.options_disp.disp_var
                if options.options_disp.disp_legend_var
                    vis_mn = 'on';
                    if options.options_disp.disp_legend_compact
                        legend_str{end+1} = ['$\pm$ $2$ Std. ', '$(', num2str(rad2deg(theta(m)), 3), '^\circ,' num2str(rad2deg(phi(m)), 3), '^\circ)$'];
                    else
                        legend_str{end+1} = ['Posterior Mean $\pm$ $2$ Std. ', '$T_{60}(\omega,', num2str(rad2deg(theta(m)), 3), '^\circ,' num2str(rad2deg(phi(m)), 3), '^\circ)$'];
                    end
                else %options.options_disp.disp_legend_var is false, display one entry only
                    if m == 1
                        legend_str{end+1} = ['Posterior Mean $\pm$ $2$ Std.'];
                        vis_mn = 'on';
                    else
                        vis_mn = 'off';
                    end
                end
                fill([freq; flipud(freq)]', [exp( (log(disp_mu_T60(:, m)) + 2 * sqrt(disp_var_log_T60(:, m)) ));  flipud(exp( (log(disp_mu_T60(:, m)) - 2 * sqrt(disp_var_log_T60(:, m)) ))) ]',  ...
                    [0.5 0.5 0.5], 'FaceAlpha', options.options_disp.disp_var_transparency, 'EdgeAlpha', 0, 'HandleVisibility', vis_mn); hold on;
            end
        end
    end

    %Plot sample evaluations
    if options.options_disp.disp_sample_eval
        for n = 1:num_evals
            for m = 1:options.options_disp.disp_sample_stride:N_E
                if options.options_disp.disp_legend_samples
                    vis_mn = 'on';
                    if options.options_disp.disp_legend_compact
                        if num_evals > 1
                            legend_str{end+1} = ['S = ', num2str(n), ': $(', num2str(rad2deg(theta(m)), 3), '^\circ,' num2str(rad2deg(phi(m)), 3), '^\circ)$'];
                        else
                            legend_str{end+1} = ['$(', num2str(rad2deg(theta(m)), 3), '^\circ,' num2str(rad2deg(phi(m)), 3), '^\circ)$'];
                        end    
                    else
                        if num_evals > 1
                            legend_str{end+1} = ['Sample ', num2str(n), ': $T_{60}(\omega,', num2str(rad2deg(theta(m)), 3), '^\circ,' num2str(rad2deg(phi(m)), 3), '^\circ)$'];
                        else
                            legend_str{end+1} = ['$T_{60}(\omega,', num2str(rad2deg(theta(m)), 3), '^\circ,' num2str(rad2deg(phi(m)), 3), '^\circ)$'];
                        end                
                    end
                else %disp_legend_samples is false, display one entry only
                    if n == 1 && m == 1
                        legend_str{end+1} = ['Sample T60'];
                        vis_mn = 'on';
                    else
                        vis_mn = 'off';
                    end
                end
                semilogx(freq, exp(log_T60(:, m, n)), 'linewidth', 2, 'HandleVisibility', vis_mn); hold on;
            end
        end
    end

    %Plot observations
    if has_obs
        semilogx(freq_obs, obs.T60, 'r*', 'MarkerSize', options.options_disp.disp_marker_size, 'linewidth', 2.5); hold on;
        legend_str{end+1} = 'Observed T60';
    end
    grid on; axis tight;

    xlabel('Frequency (Hz)', 'fontsize', fontsize);
    ylabel('T60 (Second)', 'fontsize', fontsize);
    title(['Gaussian Process ', name_mod, ' Sample T60s'], 'fontsize', fontsize + 1);
    set(gca, 'fontsize', fontsize - 1);

    h_lg = legend(legend_str, 'location', options.options_disp.disp_legend_loc, ...
        'interpreter', 'latex', 'NumColumns', options.options_disp.disp_legend_num_cols, ...
        'BackgroundAlpha', options.options_disp.disp_legend_transparency);
    set(h_lg, 'fontsize', fontsize - 1);

    %Limits
    xlim(options.options_disp.disp_xlim);
    ylim(options.options_disp.disp_ylim);

end