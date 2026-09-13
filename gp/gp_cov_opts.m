function options = gp_cov_opts(options)
%Get default GP prior covariance options

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input and output

%options:               struct
%options.cov_func:      String, covariance function name {'cov_sqx_chw_ns', 'cov_sqx_chw_sqx'})
%                           'cov_sqx_chw_ns':    Squared exponential of chordal distance with non-stationary frequency (cov_sqx_chw_ns.m)
%                           'cov_sqx_chw_sqx':   Squared exponential of chordal distance x squared exponential of frequency (cov_sqx_chw_sqx.m)

arguments    
    %Default function
    options.cov_func   (1,:) char {mustBeMember(options.cov_func, {'cov_sqx_chw_ns', 'cov_sqx_chw_sqx'}) } = 'cov_sqx_chw_ns';
    
    %Default parameters
    options.cov_sigma  (1,1) double {mustBePositive}    = sqrt(2)/2;
    options.cov_ell    (1,1) double {mustBePositive}    = 343;
    options.cov_gamma  (1,1) double  = 1;
    options.cov_jit    (1,1) double {mustBeNonnegative} = 1e-20;

    options.cov_ell_c    (1,1) double {mustBePositive}    = 1;
    options.cov_ell_f    (1,1) double {mustBePositive}    = 2;

    %Default limits [min, max]
    options.cov_sigma_lim (1,2) double {mustBePositive} = [sqrt(eps), inf];
    options.cov_ell_lim (1,2) double {mustBePositive} = [eps, inf];
    options.cov_gamma_lim (1,2) double {mustBeNonnegative}  = [0, inf];
    options.cov_jit_lim (1,2) double {mustBeNonnegative} = [0, inf];

    options.cov_ell_c_lim  (1,2) double {mustBePositive} = [eps, inf];
    options.cov_ell_f_lim  (1,2) double {mustBePositive} = [eps, inf];

end

if strcmp(options.cov_func, 'cov_sqx_chw_ns')

    assert(options.cov_sigma >= min(options.cov_sigma_lim) && options.cov_sigma <= max(options.cov_sigma_lim), 'options.cov_sigma out of bounds')
    assert(options.cov_ell >= min(options.cov_ell_lim) && options.cov_ell <= max(options.cov_ell_lim), 'options.cov_ell out of bounds')
    assert(options.cov_gamma >= min(options.cov_gamma_lim) && options.cov_gamma <= max(options.cov_gamma_lim), 'options.cov_gamma out of bounds')
    assert(options.cov_jit >= min(options.cov_jit_lim) && options.cov_jit <= max(options.cov_jit_lim), 'options.cov_jit out of bounds')

elseif strcmp(options.cov_func, 'cov_sqx_chw_sqx')

    assert(options.cov_sigma >= min(options.cov_sigma_lim) && options.cov_sigma <= max(options.cov_sigma_lim), 'options.cov_sigma out of bounds')
    assert(options.cov_ell_c >= min(options.cov_ell_c_lim) && options.cov_ell_c <= max(options.cov_ell_c_lim), 'options.cov_ell_c out of bounds')
    assert(options.cov_ell_f >= min(options.cov_ell_f_lim) && options.cov_ell_f <= max(options.cov_ell_f_lim), 'options.cov_ell_f out of bounds')
    assert(options.cov_jit >= min(options.cov_jit_lim) && options.cov_jit <= max(options.cov_jit_lim), 'options.cov_jit out of bounds')

else
    error('Unsupported options.cov_func');
end
