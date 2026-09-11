function [pass, pass_list] = tst_gp_check_grads(cov_name, options)
%Numerically check analytic gradients:  dlog(det(K)) / dparam of covariance function
%and gram matrix K

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%cov_name:          String, covariance function name {'cov_sqx_chw_ns', 'cov_sqx_chw_sqx', 'cov_sqx_chw'}

%options:           Struct
%options.rseed:     Random seed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%pass:          Logical, if true, pass all test cases 
%pass_list:     [1 x *]  Logical list for each test case

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Test gradients for various covariance functions

%pass = tst_gp_check_grads('cov_sqx_chw_ns')
%pass = tst_gp_check_grads('cov_sqx_chw_sqx')
%pass = tst_gp_check_grads('cov_sqx_chw')

arguments
    cov_name (1,:) char {mustBeMember(cov_name, {'cov_sqx_chw_ns', 'cov_sqx_chw_sqx', 'cov_sqx_chw'})} = 'cov_sqx_chw_ns';

    options.rseed (1,1) double {mustBePositive, mustBeInteger} = 213;
end

rng(options.rseed);

if strcmp(cov_name, 'cov_sqx_chw_ns')

    N = 10;
    X = [1 ./ (1:N)', rand(N, 1) * pi, rand(N, 1) * 2 * pi];
    sigma = 0.4;
    ell = 5.3;
    gamma = 1.2;

    [grad_check_sigma] = checkGradients(@(x) func_grad_test_cov_sqx_chw_ns_sigma(x, X, ell, gamma), sigma, 'Display', 'on')
    [grad_check_ell] = checkGradients(@(x) func_grad_test_cov_sqx_chw_ns_ell(x, X, sigma, gamma), ell , 'Display', 'on')
    [grad_check_gamma] = checkGradients(@(x) func_grad_test_cov_sqx_chw_ns_gamma(x, X, ell, sigma), gamma , 'Display', 'on')
    
    pass_list = [grad_check_sigma, grad_check_ell, grad_check_gamma];

elseif strcmp(cov_name, 'cov_sqx_chw_sqx')

    N = 10;
    X = [logspace(log10(50), log10(24000), N)', rand(N, 1) * pi, rand(N, 1) * 2 * pi];
    sigma = 0.4;
    ell_c = 5.3;
    ell_f = 2.3;

    [grad_check_sigma] = checkGradients(@(x) func_grad_test_cov_sqx_chw_sqx_sigma(x, X, ell_c, ell_f), sigma, 'Display', 'on')
    [grad_check_ell_c] = checkGradients(@(x) func_grad_test_cov_sqx_chw_sqx_ell_c(x, X, sigma, ell_f), ell_c , 'Display', 'on')
    [grad_check_ell_f] = checkGradients(@(x) func_grad_test_cov_sqx_chw_sqx_ell_f(x, X, sigma, ell_c), ell_f , 'Display', 'on')
    
    pass_list = [grad_check_sigma, grad_check_ell_c, grad_check_ell_f];

elseif strcmp(cov_name, 'cov_sqx_chw')

    N = 10;
    X = [rand(N, 1) * pi, rand(N, 1) * 2 * pi];
    sigma = 0.4;
    ell = 5.3;

    [grad_check_sigma] = checkGradients(@(x) func_grad_test_cov_sqx_chw_sigma(x, X, ell), sigma, 'Display', 'on')
    [grad_check_ell] = checkGradients(@(x) func_grad_test_cov_sqx_chw_ell(x, X, sigma), ell , 'Display', 'on')
    
    pass_list = [grad_check_sigma, grad_check_ell];

else
    error('Unsupported cov_name');
end

pass = all(pass_list);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Covariance function: cov_sqx_chw_ns.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fval, grad] = func_grad_test_cov_sqx_chw_ns_sigma(sigma, X, ell, gamma)
%fval = log(deg(K + eye(size(K, 1)))

[K, dK_dsigma] = cov_sqx_chw_ns(X, X, sigma, ell, gamma);
K_regu = K + eye(size(K, 1));
fval = log(det(K_regu));
grad = trace(K_regu \ dK_dsigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fval, grad] = func_grad_test_cov_sqx_chw_ns_ell(ell, X, sigma, gamma)
%fval = log(deg(K + eye(size(K, 1)))

[K, ~, dK_dell] = cov_sqx_chw_ns(X, X, sigma, ell, gamma);
K_regu = K + eye(size(K, 1));
fval = log(det(K_regu));
grad = trace(K_regu \ dK_dell);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fval, grad] = func_grad_test_cov_sqx_chw_ns_gamma(gamma, X, ell, sigma)
%fval = log(deg(K + eye(size(K, 1)))

[K, ~, ~, dK_dgamma] = cov_sqx_chw_ns(X, X, sigma, ell, gamma);
K_regu = K + eye(size(K, 1));
fval = log(det(K_regu));
grad = trace(K_regu \ dK_dgamma);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Covariance function: cov_sqx_chw_sqx.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fval, grad] = func_grad_test_cov_sqx_chw_sqx_sigma(sigma, X, ell_c, ell_f)
%fval = log(deg(K + eye(size(K, 1)))

[K, dK_dsigma] = cov_sqx_chw_sqx(X, X, sigma, ell_c, ell_f);
K_regu = K + eye(size(K, 1));
fval = log(det(K_regu));
grad = trace(K_regu \ dK_dsigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fval, grad] = func_grad_test_cov_sqx_chw_sqx_ell_c(ell_c, X, sigma, ell_f)
%fval = log(deg(K + eye(size(K, 1)))

[K, ~, dK_dell_c] = cov_sqx_chw_sqx(X, X, sigma, ell_c, ell_f);
K_regu = K + eye(size(K, 1));
fval = log(det(K_regu));
grad = trace(K_regu \ dK_dell_c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fval, grad] = func_grad_test_cov_sqx_chw_sqx_ell_f(ell_f, X, sigma, ell_c)
%fval = log(deg(K + eye(size(K, 1)))

[K, ~, ~, dK_dell_f] = cov_sqx_chw_sqx(X, X, sigma, ell_c, ell_f);
K_regu = K + eye(size(K, 1));
fval = log(det(K_regu));
grad = trace(K_regu \ dK_dell_f);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Covariance function: cov_sqx_chw.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fval, grad] = func_grad_test_cov_sqx_chw_sigma(sigma, X, ell)
%fval = log(deg(K + eye(size(K, 1)))

[K, dK_dsigma] = cov_sqx_chw(X, X, sigma, ell);
K_regu = K + eye(size(K, 1));
fval = log(det(K_regu));
grad = trace(K_regu \ dK_dsigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fval, grad] = func_grad_test_cov_sqx_chw_ell(ell, X, sigma)
%fval = log(deg(K + eye(size(K, 1)))

[K, ~, dK_dell] = cov_sqx_chw(X, X, sigma, ell);
K_regu = K + eye(size(K, 1));
fval = log(det(K_regu));
grad = trace(K_regu \ dK_dell);