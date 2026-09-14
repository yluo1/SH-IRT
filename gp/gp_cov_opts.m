function options = gp_cov_opts(options)
%Get default GP prior covariance options w.r.t. angular frequency x spherical coordinates

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input and output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments    
    %Default covariance function
    options.cov_func   (1,:) char {mustBeMember(options.cov_func, {'cov_sqx_chw_ns', 'cov_sqx_chw_sqx'}) } = 'cov_sqx_chw_ns'; %Prior covariance function name
    %   'cov_sqx_chw_ns':    Squared exponential of chordal distance with non-stationary frequency (cov_sqx_chw_ns.m)
    %   'cov_sqx_chw_sqx':   Squared exponential of chordal distance x squared exponential of frequency (cov_sqx_chw_sqx.m)
        
    %Default parameters
    options.cov_sigma  (1,1) double {mustBePositive}    = sqrt(2)/2; % Covariance scaling hyper-parameter
    options.cov_ell    (1,1) double {mustBePositive}    = 343;       % Wavelength lambda scaling (velocity) hyperparameter where lambda = ell / f^gamma
    options.cov_gamma  (1,1) double                     = 1;         % Frequency power hyper-parameter 
    options.cov_jit    (1,1) double {mustBeNonnegative} = 1e-20;     % Diagonal loading jitter: K + cov_jit * I 

    options.cov_ell_c    (1,1) double {mustBePositive}    = 1;       % Chordal distance scaling hyper-parameter
    options.cov_ell_f    (1,1) double {mustBePositive}    = 2;       % log-frequency hyper-parameter

    %Default limits [min, max]
    options.cov_sigma_lim (1,2) double {mustBePositive}     = [sqrt(eps), inf];
    options.cov_ell_lim (1,2) double {mustBePositive}       = [eps, inf];
    options.cov_gamma_lim (1,2) double {mustBeNonnegative}  = [0, inf];
    options.cov_jit_lim (1,2) double {mustBeNonnegative}    = [0, inf];

    options.cov_ell_c_lim  (1,2) double {mustBePositive}    = [eps, inf];
    options.cov_ell_f_lim  (1,2) double {mustBePositive}    = [eps, inf];

    %Misc.
    options.enable_disp (1,1) logical = false;                       % Logical, if true, display plot in gp_cov.m
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
