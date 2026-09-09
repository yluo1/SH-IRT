function tst_gp_check_grads
%Numerically check analytic gradients:  dlog(det(K)) / dparam

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(213);

N = 10;
X = [1 ./ (1:N)', rand(N, 1) * pi, rand(N, 1) * 2 * pi];
sigma = 0.4;
ell = 5.3;
gamma = 1.2;

[grad_check_sigma, err_sigma] = checkGradients(@(x) func_grad_test_sigma(x, X, ell, gamma), sigma, 'Display', 'on')
[grad_check_ell, err_sigma] = checkGradients(@(x) func_grad_test_ell(x, X, sigma, gamma), ell , 'Display', 'on')
[grad_check_gamma, err_gamma] = checkGradients(@(x) func_grad_test_gamma(x, X, ell, sigma), gamma , 'Display', 'on')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fval, grad] = func_grad_test_sigma(sigma, X, ell, gamma)
%fval = log(deg(K + eye(size(K, 1)))

[K, dK_dsigma] = cov_sqx_chw_ns(X, X, sigma, ell, gamma);
K_regu = K + eye(size(K, 1));
fval = log(det(K_regu));
grad = trace(K_regu \ dK_dsigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fval, grad] = func_grad_test_ell(ell, X, sigma, gamma)
%fval = log(deg(K + eye(size(K, 1)))

[K, ~, dK_dell] = cov_sqx_chw_ns(X, X, sigma, ell, gamma);
K_regu = K + eye(size(K, 1));
fval = log(det(K_regu));
grad = trace(K_regu \ dK_dell);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fval, grad] = func_grad_test_gamma(gamma, X, ell, sigma)
%fval = log(deg(K + eye(size(K, 1)))

[K, ~, ~, dK_dgamma] = cov_sqx_chw_ns(X, X, sigma, ell, gamma);
K_regu = K + eye(size(K, 1));
fval = log(det(K_regu));
grad = trace(K_regu \ dK_dgamma);
