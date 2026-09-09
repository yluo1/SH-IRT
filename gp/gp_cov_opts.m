function options = gp_cov_opts(options)
%Get default GP prior covariance options

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input and output
arguments    
    %Default function
    options.cov_func   (1,:) char {mustBeMember(options.cov_func, {'SqChordalDistFreq'}) } = 'SqChordalDistFreq';
    
    %Default parameters
    options.cov_sigma  (1,1) double {mustBePositive}    = sqrt(2)/2;
    options.cov_ell    (1,1) double {mustBePositive}    = 343;
    options.cov_gamma  (1,1) double  = 1;
    options.cov_jit    (1,1) double {mustBeNonnegative} = 1e-20;

    %Default limits [min, max]
    options.cov_sigma_lim (1,2) double {mustBePositive} = [sqrt(eps), inf];
    options.cov_ell_lim (1,2) double {mustBePositive} = [eps, inf];
    options.cov_gamma_lim (1,2) double {mustBeNonnegative}  = [0, inf];
    options.cov_jit_lim (1,2) double {mustBeNonnegative} = [0, inf];

end

assert(options.cov_sigma >= min(options.cov_sigma_lim) && options.cov_sigma <= max(options.cov_sigma_lim), 'options.cov_sigma out of bounds')
assert(options.cov_ell >= min(options.cov_ell_lim) && options.cov_ell <= max(options.cov_ell_lim), 'options.cov_ell out of bounds')
assert(options.cov_gamma >= min(options.cov_gamma_lim) && options.cov_gamma <= max(options.cov_gamma_lim), 'options.cov_gamma out of bounds')
assert(options.cov_jit >= min(options.cov_jit_lim) && options.cov_jit <= max(options.cov_jit_lim), 'options.cov_jit out of bounds')

