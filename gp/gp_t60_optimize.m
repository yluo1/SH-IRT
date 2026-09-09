function [options_cov, lmh] = gp_t60_optimize(obs, options)
%Optimize gaussian process T60(omega, theta, phi) hyperparameters via
%maximum log-marginal likelihood

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%obs:                   Struct, observed T60(omega, theta, phi), see gp_obs_opts.m

%options:               Struct, 
%options.options_mu:    Struct, prior mean options, see gp_mu_opts.m
%options.options_cov:   Struct, prior covariance options, see gp_cov_opts.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%options_cov:           Struct, posterior covariance options, see gp_cov_opts.m
%lmh:                   Log-marginal likelihood

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:
% rng(1278);
% N_S = 200;
% [theta, phi] = sh_fib(N_S);
% omega = 1000 * 2 * pi * ones(N_S, 1);
% T60 = rand(N_S, 1) * 0.1 + 0.5;
% log_noise_std = ones(N_S ,1) * 0.01; %1 percent
% obs = gp_obs_opts('omega', omega, 'theta', theta, 'phi', phi, 'T60', T60, 'log_noise_std', log_noise_std);

%gp_t60_optimize(obs);

arguments

    %Observations
    obs = gp_obs_opts('omega', 1000*[2*pi 2*pi]', 'theta', deg2rad([90, 90])', 'phi', deg2rad([0, 30])', 'T60', [1.1 0.9]', 'log_noise_std', [1e-5 1e-5]');

    %Priors
    options.options_mu  = gp_mu_opts;
    options.options_cov = gp_cov_opts;

    options.options_fmincon = optimoptions("fmincon", SpecifyObjectiveGradient=true, Display="iter", checkGradients=false, ...
        MaxIterations=1000);

    %Display
    options.enable_disp (1,1) logical = false;
    options.options_disp = gp_disp_opts;

end

options_cov = options.options_cov;

%Setup input set {X, y}
freq_obs = max(1, obs.omega / (2 * pi) ); %[N_S x 1]
lambda_obs = 1 ./ freq_obs;
X = [lambda_obs(:), obs.theta(:), obs.phi(:)];
y = log(obs.T60) - log(gp_mu(lambda_obs(:), options.options_mu));

%Get number of parameters
[~, ~, dK_name_list] = gp_cov(X, X, options_cov);
N_params = numel(dK_name_list);

%Setup fmincon
param_list_0 = zeros(N_params, 1);
lb = -inf(N_params, 1);
ub = inf(N_params, 1);
for n = 1:N_params
    param_list_0(n) = getfield(options_cov, dK_name_list{n});
    lim_n = getfield(options_cov, [dK_name_list{n}, '_lim']);
    lb(n) = min(lim_n);
    ub(n) = max(lim_n);
end

%Optimize
[param_list, fval] = fmincon(@(x) neg_log_marginal_likelihood(x, X, y, obs.log_noise_std, options_cov, dK_name_list), ... 
    param_list_0, [], [], [], [], lb, ub, [], options.options_fmincon);

lmh = -fval;

%Write to output struct
disp(['log-marginal likelihood: ', num2str(lmh)]);
for n = 1:N_params
    options_cov = setfield(options_cov, dK_name_list{n}, param_list(n));
    disp([dK_name_list{n}, ': ', num2str(param_list(n))]);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fval, grad] = neg_log_marginal_likelihood(param_list, X, y, log_noise_std, options_cov, dK_name_list)

N_params = numel(dK_name_list);
for n = 1:N_params
    options_cov = setfield(options_cov, dK_name_list{n}, param_list(n));
end

[K, dK_mat_list] = gp_cov(X, X, options_cov);
K = K + diag(log_noise_std.^2);

%Compute negative log-marginal likelihood
a = K \ y;
fval = 1/2 * (y'*a + log(det(K)) + numel(y) * log(2*pi) );

%Compute partial derivative w.r.t. each variable
grad = zeros(N_params, 1);
for n = 1:N_params
    grad(n) = 1/2 * trace((K \ dK_mat_list{n}) - a*a' * dK_mat_list{n} );
end

