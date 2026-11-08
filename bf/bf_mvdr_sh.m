function [W, fval, SINR, WNG, h_figs] = bf_mvdr_sh(C_arr, is_real, freq, theta, phi, kappa, mode, options)
% Minimum variance distortionless response beamformer designs (bf_mvdr.m)
% with steering vectors and spatial covariances R computed from spherical harmonic expansions C

% mode = 'uniform'      Compute R over uniform density
% mode = 'density'      Compute R over density options.R_pdf
% mode = 'direction'    Compute R over steering-directions options.R_theta_phi 
% mode = 'custom'       Specify R in options.R_custom

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C_arr:     [(L_C + 1)^2 X N x M]   SH expansion of array, L_C max-order, N channels, M frequencies
%is_real:   Logical, if true, evaluate real SH

%freq:      [1 x M]   Frequency in Hz

%theta:     Steering co-latitude (0, pi)
%phi:       Steering azimuth (0, 2 * pi)

%kappa:     Scalar, condition number constraint >= 1, inf = unconstrained
%mode:      String, spatial covariance compute method {'uniform', 'density', 'direction', 'custom'}

%options:               struct
%options.mode_kappa:    String, condition number constraint method {'diag', 'thres', 'exp', 'sse', 'ase', 'minimax'} (see bf_mvdr.m)

%options.R_pdf:         [(L_R + 1)^2 X 1] or [(L_R + 1)^2 X M]  SH density expansion of L_R max-order over M frequencies.
%                                                               If second dimension is 1, duplicate expansions over M frequencies
%options.R_theta_phi:   [N_theta_phi x 2]  Steering directions (theta, phi) to compute R for mode = 'direction'
%options.R_custom:      [N x N x M] Custom covariance matrix over frequency for mode = 'custom'

%options.center:    [] or [N x 3]   Center-of-expansion coordinates xyz per channel 
%options.wavnum:    [] or [1 x M]   Angular wavenumber kappa (e.g. kappa = 2 * pi * hz / c for propagation velocity c)

%options.wt:        [] or [N x M]   Frequency weighting per channel / frequency. If empty, treat channels as equal-weighted

%options.enable_disp:   Logical, if true, display beampattern and statistics
%options.disp_hplane:   Logical, if true, display horizontal plane beam pattern
%options.disp_mplane:   Logical, if true, display median plane beam pattern
%options.disp_fplane:   Logical, if true, display frontal plane beam pattern
%options.disp_stats:    Logical, if true, display stats
%options.disp_dB_lim:   [1 x 2]  Display dB limits [min, max]
%options.disp_freq:     [1 x M_f]   Display beampatterns at select frequencies if non-empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%W:         [N x M]     Beamformer weights (conjugated)
%fval       [1 x M]     Objective:           fval(m) = real(W(:,m)' * R_kappa(:,:,m) * W(:,m))
%SINR:      [1 x M]     Signal-to-interference-plus-noise ratio:   SINR(m) = |W(:,m)' * v_m|^2 / real(W(:,m)' * R_m * W(:,m))
%WNG:       [1 x M]     White-noise gain:    WNG(m) = W(:,m)' * W(:,m)
%h_figs:    [1 x *]     Cell array figure handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Design MVDR beamformer over circular omni-directivity transducer array

% Fs = 16000;
% M = 256;
% N = 3;
% radius = 0.1;
% delay_dist = 1;
% max_odr = 36;
% mode = 'uni_circ';
% freq = logspace(log10(20), log10(Fs), M + 1); freq = freq(1:end-1);

% theta = pi/2;
% phi = 0;
% kappa = inf;
% mode_mvdr = 'uniform';

% C_arr = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, 'enable_disp', true, 'disp_dB_lim', [-3, 3]);
% bf_mvdr_sh(C_arr, false, freq, theta, phi, kappa, mode_mvdr, 'enable_disp', true, 'disp_dB_lim', [-48, 14]);

% [C_arr_noncenter, center] = sh_enc_arr(0, mode, N, radius, delay_dist, freq, 'noncentered', true, 'enable_disp', true, 'disp_dB_lim', [-3, 3]);
% bf_mvdr_sh(C_arr_noncenter, false, freq, theta, phi, kappa, mode_mvdr, 'center', center, 'wavnum', 2 * pi * freq / 343, 'enable_disp', true, 'disp_dB_lim', [-48, 14]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Design MVDR beamformer over circular piston directivity transducer array

% Fs = 24000;
% M = 256;
% N = 3;
% radius = 0.0254 * 3; % 3 inch
% delay_dist = 1;
% max_odr = 30;
% mode = 'spis_circ';
% freq = logspace(log10(20), log10(Fs), M + 1); freq = freq(1:end-1);
% spis_arc_len = 0.0254;

% theta = pi/2;
% phi = 0;
% kappa = inf;
% mode_mvdr = 'uniform';

% C_arr = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, 'spis_arc_len', spis_arc_len, 'enable_disp', true, 'disp_dB_lim', [-36, 12]);
% bf_mvdr_sh(C_arr, false, freq, theta, phi, kappa, mode_mvdr, 'enable_disp', true, 'disp_dB_lim', [-48, 14]);


arguments
    C_arr (:,:,:) double {coder.mustBeComplex} = complex(0)
    is_real (1,1) logical = false;
            
    freq (1,:) double {mustBeNonnegative} = 0;

    theta (1,1) double = pi/2;
    phi (1,1) double = 0;

    kappa (1,1) double {mustBePositive} = inf;
    mode (1,:) char {mustBeMember(mode, {'uniform', 'density', 'direction', 'custom'} )} = 'uniform';

    options.mode_kappa (1,:) char {mustBeMember(options.mode_kappa, {'diag', 'thres', 'exp', 'sse', 'ase', 'minimax'} )} = 'diag';

    options.R_pdf (:,:) double {coder.mustBeComplex} = sh_pdf_preset('FwdCenter', 'Uniform', 0);
    options.R_theta_phi (:,2) double = [];
    options.R_custom (:,:,:) double  {coder.mustBeComplex} = [];

    options.center (:,3) double = [];
    options.wavnum (1,:) double = [];

    options.wt (:,:) double {mustBeNonnegative}  = [];

    options.enable_disp (1,1) logical = false;
    options.disp_hplane (1,1) logical = true;
    options.disp_mplane (1,1) logical = false;
    options.disp_fplane (1,1) logical = false;
    options.disp_stats (1,1) logical = true;
    options.disp_dB_lim (1,2) double = [-inf, inf];
    options.disp_freq (1,:) double {mustBeNonnegative} = [];

end

[N_C, N, M] = size(C_arr);
L_C = sqrt(N_C) - 1;
assert( floor(L_C) == L_C, 'Invalid size C');

if ~isempty(options.center)
    assert( size(options.center, 1) == N, 'Size options.center mismatch N channels');
    assert( ~isempty(options.wavnum), 'Requires non-empty options.wavnum' );
end
if ~isempty(options.wavnum)
    assert( numel(options.wavnum) == M, 'Size options.wavnum mismatch M frequencies');
end

if ~isempty(options.wt)
    assert( isequal(size(options.wt), [N, M]), 'options.wt mismatch [N, M]');
end

% Compute spatial covariance matrix
if strcmp(mode, 'uniform')

    R = sh_cov(C_arr, sh_pdf_preset('FwdCenter', 'Uniform', 0), is_real, 'SHT_center', options.center, 'SHT_wavnum', options.wavnum, ...
        'SHT_oversample', 1, 'SHT_maxsample', 4096); % [N x N x M] Spatial covariance integrated over uniform density

elseif strcmp(mode, 'density')

    assert(~isempty(options.R_pdf), 'options.R_pdf must be non-empty');

    R = sh_cov(C_arr, options.R_pdf, is_real, 'SHT_center', options.center, 'SHT_wavnum', options.wavnum, ...
        'SHT_oversample', 1, 'SHT_maxsample', 4096); % [N x N x M] Spatial covariance integrated over density function R_pdf

elseif strcmp(mode, 'direction')

    assert(~isempty(options.R_theta_phi), 'options.R_theta_phi must be non-empty');
    R = complex(zeros([N, N, M]));  % [N x N x M] Spatial covariance along directions options.R_theta_phi

    for m = 1:M

        if isempty(options.center)
            center_wavnum = [];
        else
            center_wavnum = [options.center, options.wavnum(m) * ones(N,1)];
        end
        V = sh_dec(C_arr(:, :, m), options.R_theta_phi(:, 1), options.R_theta_phi(:, 2), is_real, 'eval', center_wavnum).'; % [N x N_theta_phi]

        R(:,:,m) = V*V'; 
    end
elseif strcmp(mode, 'custom')

    assert(isequal(size(options.R_custom), [N, N, M]), 'options.R_custom size mismatch [N x N x M]');

    R = options.R_custom;

else
    error('Unknown mode');
end

% Compute beamformer weights
W = complex(zeros(N, M));
fval = zeros(1,M);
SINR = zeros(1,M);
WNG = zeros(1,M);

for m = 1:M

    % Channel / Frequency weighting
    if isempty(options.wt) 
        R_m = R(:,:,m);
    else 
        R_m = diag(options.wt(:,m)) * R(:,:,m) * diag(options.wt(:,m)); 
    end

    % Steering vector
    if isempty(options.center)
        center_wavnum = [];
    else
        center_wavnum = [options.center, options.wavnum(m) * ones(N,1)];
    end
    v_m = sh_dec(C_arr(:, :, m), theta, phi, is_real, 'eval', center_wavnum).';  % [N x 1] Steering vector

    % MVDR solver
    [W(:, m), fval(m), SINR(m), WNG(m)] = bf_mvdr(v_m, R_m, kappa, options.mode_kappa);

end

% Display
h_figs = [];
if options.enable_disp && coder.target("MATLAB")

    varargin_bf_plt_sh = {'center', options.center, 'wavnum', options.wavnum, 'dB_lim', options.disp_dB_lim, 'title_name', 'MVDR', ...
        'disp_freq', options.disp_freq};

    if options.disp_hplane
        h_fig_bf = bf_plt_sh(C_arr, false, W, freq, 'horizontal', varargin_bf_plt_sh{:});
        h_figs{end+1} = h_fig_bf;
    end

    if options.disp_mplane
        h_fig_bf = bf_plt_sh(C_arr, false, W, freq, 'median', varargin_bf_plt_sh{:});
        h_figs{end+1} = h_fig_bf;
    end

    if options.disp_fplane
        h_fig_bf = bf_plt_sh(C_arr, false, W, freq, 'frontal', varargin_bf_plt_sh{:});
        h_figs{end+1} = h_fig_bf;
    end

    if options.disp_stats
        title_name = ['mode = ', mode, ', \kappa = (', options.mode_kappa, ', ', num2str(kappa), '), direction (\theta = ', num2str(rad2deg(theta)), ', \phi = ', num2str(rad2deg(phi)), ')'];
        h_fig_stats = bf_plt_sh(C_arr, false, W, freq, 'stats', varargin_bf_plt_sh{:}, 'title_name', title_name);
        h_figs{end+1} = h_fig_stats;
    end
    
end