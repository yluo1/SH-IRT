function [W, h_figs] = bf_grq_sh(C_arr, is_real, freq, tau_tgt, mode_sol, mode_A, mode_R, options)
% Generalized Rayleigh quotient (GRQ) beamformer designs for spherical harmonic expansions 
% constrain the solution to smooth and bounded generalized directivity index functions
% GRQ(A, R, W) = tau \in [min gev, max gev] w.r.t. frequency.

% mode_sol = 'spec'     Minimum weighted spectral constant directivity index (bf_grq_spec.m)

% mode_sol = 'proj'     Minimum norm projection onto constant directivity constraint (bf_grq_proj.m)

% mode_sol = 'sensitivity'      Maximum sensitivity constant directivity (bf_grq_mscd.m)
% mode_sol = 'sensitivity-nowt' Maximum sensitivity constant non-weighted regularized channel directivity (bf_grq_mscd.m)

% mode_sol = 'efficiency'      Maximum efficiency constant directivity  (bf_grq_mecd.m)
% mode_sol = 'efficiency-nowt' Maximum efficiency constant non-weighted channel directivity  (bf_grq_mecd.m)

% mode_{A,R} = 'uniform'      Compute {A,R} over uniform density
% mode_{A,R} = 'density'      Compute {A,R} over SH density expansion in options.{A,R}_pdf
% mode_{A,R} = 'direction'    Compute {A,R} over spherical coordinates in options.{A,R}_theta_phi 
% mode_{A,R} = 'custom'       Specify {A,R} in options.{A,R}_custom

%Paper References:
%Luo, Y. (2021). Spherical harmonic covariance and magnitude function encodings for beamformer design. EURASIP Journal on Audio, Speech, and Music Processing, 2021(1), 41. 
%Luo, Y. (2024). Constant Directivity Loudspeaker Beamforming. In 2024 32nd European Signal Processing Conference (EUSIPCO) (pp. 246-250). IEEE.

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C_arr:         [(L_C + 1)^2 x N x M]   SH coefficients (max order P of N number of channels expansions, M frequencies)

%is_real:       Logical, if true, evaluate real SH

%freq:          [1 x M]         Frequency (Hz)
%tau_tgt:       Scalar          Constant target                           
%               [1 x M]         Target per frequency

%mode_sol:      String, solver method   {'spec', 'proj', 'sensitivity', 'sensitivity-nowt', 'efficiency', 'efficiency-nowt'}
%mode_A:        String, spatial covariance compute method for covariance matrix A {'uniform', 'density', 'direction', 'custom'}
%mode_A:        String, spatial covariance compute method for covariance matrix R {'uniform', 'density', 'direction', 'custom'}

%options:       Struct

%options.theta_steer:   Scalar, co-latitude (0, pi) steering direction for mode_sol = {'proj', 'sensitivity', 'sensitivity-nowt'}
%options.phi_steer:     Scalar, azimuth (0, 2 * pi) steering direction for mode_sol = {'proj', 'sensitivity', 'sensitivity-nowt'}

%options.center:    [] or [N x 3]   Center-of-expansion coordinates xyz per channel 
%options.wavnum:    [] or [1 x M]   Angular wavenumber kappa (e.g. kappa = 2 * pi * hz / c for propagation velocity c)

%options.A_pdf:      [(L_A + 1)^2 x 1]       Density SH coefficients for numerator (accept)
%options.R_pdf:      [(L_R + 1)^2 x 1]       Density SH coefficients for denominator (reject)

%options.{A,R}_theta_phi:   [N_{A,R}_theta_phi x 2]  Steering directions (theta, phi) to compute {A,R} for mode_{A,R} = 'direction'
%options.{A_R}_custom:      [N x N x M] Custom covariance matrix over frequency for mode_{A,R} = 'custom'

%options.wt:            [] or [N x M]   Frequency weighting per channel / frequency. If empty, treat channels as equal-weighted

%options.tau_lambda:    Scalar, regularization term for smoothing bounded GRQ target tau
%options.tau_margin:    Scalar, margin size to shrink tau limits in bounded GRQ target tau 

%options.mode_spec:     String, averaging method {'max', 'sq', 'wsq'} for mode_sol = 'spec'

%options.enable_disp:       Logical, if true, display beampattern and statistics
%options.disp_hplane:       Logical, if true, display horizontal plane beam pattern
%options.disp_mplane:       Logical, if true, display median plane beam pattern
%options.disp_fplane:       Logical, if true, display frontal plane beam pattern
%options.disp_stats:        Logical, if true, display directivity index, electrical gain statistics
%options.disp_dB_lim:       [1 x 2]  Display beampattern dB limits [min, max]
%options.disp_tau_dB_lim:   [1 x 2]  Display tau index dB limits [min, max]
%options.disp_DI_dB_lim:    [1 x 2]  Display directivity index dB limits [min, max]
%options.disp_EP_dB_lim:    [1 x 2]  Display electrical power dB limits [min, max]
%options.disp_freq:         [1 x M_f]   Display beampatterns at select frequencies if non-empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%W:             [N x M]     Conjugated beamformer weights per channel over frequency
%h_figs:        [1 x *]     Cell array of figure handles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Compare GRQ beamformer designs

% Fs = 24000;
% M = 512;
% N = 4;
% radius = 0.0254 * 3; % 3 inch
% delay_dist = 1;
% max_odr = 30;
% mode = 'spis_circ';
% freq = logspace(log10(20), log10(Fs), M + 1); freq = freq(1:end-1);
% %freq = linspace(0, Fs, M + 1); freq = freq(1:end-1);
% spis_arc_len = 0.0254;
% is_real = false;
% 
% C_arr = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, 'spis_arc_len', spis_arc_len, 'enable_disp', true, 'disp_dB_lim', [-36, 12]);
% C_arr_couple = cat(2, C_arr(:,1,:), C_arr(:,2,:)+C_arr(:,4,:),  C_arr(:,3,:));

% tau_tgt = db2pow(5);
% A_pdf = sh_pdf_preset('FwdCenter', 'SqExp', 8, is_real, 'ell', 0.5);
% varargin = {'A_pdf', A_pdf, 'enable_disp', true, 'disp_dB_lim', [-48, 14], 'disp_tau_dB_lim', [-32, 12], 'disp_DI_dB_lim', [-12, 24], 'disp_EP_dB_lim', [-48, 12]};

% bf_grq_sh(C_arr, false, freq, tau_tgt, 'spec', varargin{:});
% bf_grq_sh(C_arr, false, freq, tau_tgt, 'proj', varargin{:});
% bf_grq_sh(C_arr, false, freq, tau_tgt, 'efficiency', varargin{:});
% bf_grq_sh(C_arr, false, freq, tau_tgt, 'sensitivity', varargin{:});

% bf_grq_sh(C_arr_couple, false, freq, tau_tgt, 'spec', varargin{:});
% bf_grq_sh(C_arr_couple, false, freq, tau_tgt, 'proj', varargin{:});
% bf_grq_sh(C_arr_couple, false, freq, tau_tgt, 'efficiency', varargin{:});
% bf_grq_sh(C_arr_couple, false, freq, tau_tgt, 'sensitivity', varargin{:});

% %Frequency weighted array elements
% wt = ft_freq_wt_chs(freq, [4000 8000; 100, 10000; 0 1000], [24 0; 12 6; 24 24], 'enable_disp', true, 'disp_dB_lim', [-48, 6], 'soft_knee', true, 'norm_sum', true);
% bf_grq_sh(C_arr_couple, false, freq, db2pow(inf), 'spec', varargin{:});
% bf_grq_sh(C_arr_couple, false, freq, db2pow(inf), 'efficiency', varargin{:}, 'wt', wt);
% bf_grq_sh(C_arr_couple, false, freq, db2pow(inf), 'efficiency-nowt', varargin{:}, 'wt', wt);
% bf_grq_sh(C_arr_couple, false, freq, db2pow(inf), 'sensitivity-nowt', varargin{:}, 'wt', wt);

% bf_grq_sh(C_arr_couple, false, freq, db2pow(6), 'efficiency', varargin{:}, 'wt', wt);
% bf_grq_sh(C_arr_couple, false, freq, db2pow(6), 'efficiency-nowt', varargin{:}, 'wt', wt);
% bf_grq_sh(C_arr_couple, false, freq, db2pow(6), 'sensitivity', varargin{:}, 'wt', wt);
% bf_grq_sh(C_arr_couple, false, freq, db2pow(6), 'sensitivity-nowt', varargin{:}, 'wt', wt);

arguments
    C_arr (:,:,:) double {coder.mustBeComplex} = complex(0)

    is_real (1,1) logical = false;
    freq (1,:) double {mustBeNonnegative} = 0;

    tau_tgt (1,:) {mustBeNonnegative} = db2pow(6);

    mode_sol (1,:) char {mustBeMember(mode_sol, {'spec', 'proj', 'sensitivity', 'sensitivity-nowt', 'efficiency', 'efficiency-nowt'} )} = 'sensitivity';
    mode_A   (1,:) char {mustBeMember(mode_A, {'uniform', 'density', 'direction', 'custom'} )} = 'density';
    mode_R   (1,:) char {mustBeMember(mode_R, {'uniform', 'density', 'direction', 'custom'} )} = 'uniform';

    options.theta_steer (1,1) double = pi/2;
    options.phi_steer (1,1) double = 0;

    options.center (:,3) double = [];
    options.wavnum (1,:) double = [];

    options.A_pdf (:,1) double {coder.mustBeComplex} = sh_pdf_preset('FwdCenter', 'SqExp', 8, is_real, 'ell', 0.5);
    options.A_theta_phi (:,2) double = [];
    options.A_custom (:,:,:) double  {coder.mustBeComplex} = [];

    options.R_pdf (:,1) double {coder.mustBeComplex} = [];
    options.R_theta_phi (:,2) double = [];
    options.R_custom (:,:,:) double  {coder.mustBeComplex} = [];

    options.wt (:,:) double {mustBeNonnegative}  = [];

    options.tau_lambda (1,1) double {mustBeNonnegative} = 1e4;
    options.tau_margin (1,1) double {mustBeNonnegative} = 1e-8;

    options.mode_spec (1,:) char {mustBeMember(options.mode_spec, {'max', 'sq', 'wsq'} )} = 'max';

    options.A_min_eig (1,1) double {mustBeNonnegative} = 1e-8;
    options.R_min_eig (1,1) double {mustBeNonnegative}= 1e-8;

    options.enable_disp (1,1) logical = false;
    options.disp_hplane (1,1) logical = true;
    options.disp_mplane (1,1) logical = false;
    options.disp_fplane (1,1) logical = false;
    options.disp_stats (1,1) logical = true;
    options.disp_dB_lim (1,2) double = [-inf, inf];
    options.disp_tau_dB_lim (1,2) double = [-inf, inf];
    options.disp_DI_dB_lim (1,2) double = [-inf, inf];
    options.disp_EP_dB_lim (1,2) double = [-inf, inf];
    options.disp_freq (1,:) double {mustBeNonnegative} = [];

end

[N_C, N, M] = size(C_arr);
L_C = sqrt(N_C) - 1;
assert(floor(L_C) == L_C, 'Invalid size C');

% Compute spatial covariance matrices
if strcmp(mode_A, 'uniform')

    A = sh_cov(C_arr, sh_pdf_preset('FwdCenter', 'Uniform', 0), is_real, 'SHT_center', options.center, 'SHT_wavnum', options.wavnum, ...
        'SHT_oversample', 1, 'SHT_maxsample', 4096, 'min_eig', options.A_min_eig); % [N x N x M] Spatial covariance integrated over uniform density

elseif strcmp(mode_A, 'density')

    assert(~isempty(options.A_pdf), 'options.A_pdf must be non-empty');

    A = sh_cov(C_arr, options.A_pdf, is_real, 'SHT_center', options.center, 'SHT_wavnum', options.wavnum, ...
        'SHT_oversample', 1, 'SHT_maxsample', 4096, 'min_eig', options.A_min_eig); % [N x N x M] Spatial covariance integrated over density function A_pdf

elseif strcmp(mode_A, 'direction')

    assert(~isempty(options.A_theta_phi), 'options.A_theta_phi must be non-empty');
  
    A = sh_cov(C_arr, [], is_real, 'mode', 'discrete', 'D_theta_phi', options.A_theta_phi, ...
      'SHT_center', options.center, 'SHT_wavnum', options.wavnum, ...
      'SHT_oversample', 1, 'SHT_maxsample', 4096, 'min_eig', options.A_min_eig); % [N x N x M] Spatial covariance integrated over discrete directions A_theta_phi

elseif strcmp(mode_A, 'custom')

    assert(isequal(size(options.A_custom), [N, N, M]), 'options.A_custom size mismatch [N x N x M]');

    A = options.R_custom;

else
    error('Unknown mode_A');
end

if strcmp(mode_R, 'uniform')

    R = sh_cov(C_arr, sh_pdf_preset('FwdCenter', 'Uniform', 0), is_real, 'SHT_center', options.center, 'SHT_wavnum', options.wavnum, ...
        'SHT_oversample', 1, 'SHT_maxsample', 4096, 'min_eig', options.R_min_eig); % [N x N x M] Spatial covariance integrated over uniform density

elseif strcmp(mode_R, 'density')

    assert(~isempty(options.R_pdf), 'options.R_pdf must be non-empty');

    R = sh_cov(C_arr, options.R_pdf, is_real, 'SHT_center', options.center, 'SHT_wavnum', options.wavnum, ...
        'SHT_oversample', 1, 'SHT_maxsample', 4096, 'min_eig', options.R_min_eig); % [N x N x M] Spatial covariance integrated over density function R_pdf

elseif strcmp(mode_R, 'direction')

    assert(~isempty(options.R_theta_phi), 'options.R_theta_phi must be non-empty');
  
    R = sh_cov(C_arr, [], is_real, 'mode', 'discrete', 'D_theta_phi', options.R_theta_phi, ...
      'SHT_center', options.center, 'SHT_wavnum', options.wavnum, ...
      'SHT_oversample', 1, 'SHT_maxsample', 4096, 'min_eig', options.R_min_eig); % [N x N x M] Spatial covariance integrated over discrete directions R_theta_phi

elseif strcmp(mode_R, 'custom')

    assert(isequal(size(options.R_custom), [N, N, M]), 'options.R_custom size mismatch [N x N x M]');

    R = options.R_custom;

else
    error('Unknown mode_R');
end

% Non-weighted covariances
A_nowt = A;
R_nowt = R;

% Compute frequency weighted covariances
if ~isempty(options.wt)
    assert( isequal(size(options.wt), [N, M]), 'options.wt mismatch [N, M]');
    wt = options.wt;

    for m = 1:M

        wt_mat_m = wt(:,m) * wt(:,m)';
        A(:,:,m) = A(:,:,m) .* wt_mat_m;         
        R_diag_m = diag(R(:,:,m));
        R(:,:,m) = R(:,:,m) .* wt_mat_m;
        for n = 1:N % Replace diagonal with original
            R(n,n,m) = R_diag_m(n);
        end

    end
else
    wt = [];
end

% tau target options
varargin_tau_tgt = {'lambda', options.tau_lambda, 'margin', options.tau_margin, 'enable_disp', options.enable_disp, 'disp_dB_lim', options.disp_tau_dB_lim};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute beamformer weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W = complex(zeros(N, M));

if strcmp(mode_sol, 'spec')

    % Compute beamformer target
    [tau, h_fig_tau] = bf_grq_tau_tgt(A, R, freq, tau_tgt, varargin_tau_tgt{:});

    % Compute beamformer weights
    for m = 1:M
        W(:, m) = bf_grq_spec(A(:,:,m), R(:,:,m), tau(m), options.mode_spec);
    end

elseif strcmp(mode_sol, 'proj')
    
    % Compute beamformer target
    [tau, h_fig_tau] = bf_grq_tau_tgt(A, R, freq, tau_tgt, varargin_tau_tgt{:});

    % Compute beamformer weights
    for m = 1:M

        % Steering vector
        if isempty(options.center)
            center_wavnum = [];
        else
            center_wavnum = [options.center, options.wavnum(m) * ones(N,1)];
        end
        v_m = sh_dec(C_arr(:, :, m), options.theta_steer, options.phi_steer, is_real, 'eval', center_wavnum).';  % [N x 1] Steering vector
        W(:, m) = bf_grq_proj(v_m, A(:,:,m), R(:,:,m), tau(m), 'mode', 'secular');

    end

elseif strcmp(mode_sol, 'sensitivity')

    % Compute beamformer target
    [tau, h_fig_tau] = bf_grq_tau_tgt(A, R, freq, tau_tgt, varargin_tau_tgt{:});

    % Compute beamformer weights
    for m = 1:M
        
        % Steering vector
        if isempty(options.center)
            center_wavnum = [];
        else
            center_wavnum = [options.center, options.wavnum(m) * ones(N,1)];
        end
        v_m = sh_dec(C_arr(:, :, m), options.theta_steer, options.phi_steer, is_real, 'eval', center_wavnum).';  % [N x 1] Steering vector
        W(:, m) = bf_grq_mscd(v_m, A(:,:,m), R(:,:,m), tau(m));

    end

elseif strcmp(mode_sol, 'sensitivity-nowt')

    assert(~isempty(wt), 'Requires non-empty options.wt');

    % Compute beamformer weights under change-of-variables
    % W = diag(wt) * Y with additional regularization of covariances:
    % W'*A*W / (W'*R*W) = tau   =>  Y'* A_WT * Y / (Y'* R_WT * Y),
    % A_WT = diag(wt) * A * diag(wt) + lambda_A * eye(N),
    % R_WT = diag(wt) * R * diag(wt) + lambda_R * eye(N),

    A_WT = A_nowt;
    R_WT = R_nowt;
    
    for m = 1:M
        wt_mat = wt(:,m) * wt(:,m)';
        A_WT(:,:,m) = A_WT(:,:,m) .* wt_mat;
        R_WT(:,:,m) = R_WT(:,:,m) .* wt_mat;

        if options.A_min_eig > 0 % Add diagonal loading A_min_eig to A
            d = eigs(A_WT(:, :, m), 1, 'smallestreal');
            if d < options.A_min_eig
                A_WT(:, :, m) = A_WT(:, :, m) + (options.A_min_eig - d) * eye(N);
            end
        end
        if options.R_min_eig > 0 % Add diagonal loading R_min_eig to R
            d = eigs(R_WT(:, :, m), 1, 'smallestreal');
            if d < options.R_min_eig
                R_WT(:, :, m) = R_WT(:, :, m) + (options.A_min_eig - d) * eye(N);
            end
        end
    end

    % Compute beamformer target
    [tau, h_fig_tau] = bf_grq_tau_tgt(A_WT, R_WT, freq, tau_tgt, varargin_tau_tgt{:});
    
    % Compute beamformer weights
    for m = 1:M
        
        % Steering vector
        if isempty(options.center)
            center_wavnum = [];
        else
            center_wavnum = [options.center, options.wavnum(m) * ones(N,1)];
        end
        v_m = sh_dec(C_arr(:, :, m), options.theta_steer, options.phi_steer, is_real, 'eval', center_wavnum).';  % [N x 1] Steering vector
        
        Y_m = bf_grq_mscd(v_m .* wt(:,m), A_WT(:,:,m), R_WT(:,:,m), tau(m));
        W(:, m) = Y_m .* wt(:,m); % Recover beamformer weights
        % Check constraints at solution
        % abs(W(:, m)'*A_nowt(:,:,m)*W(:,m) / ( W(:, m)'*R_nowt(:,:,m)*W(:,m) ) - tau(m))
        % W(:,m)' * w_m
        ;
    end

    % Revert A, R for display
    A = A_nowt;
    R = R_nowt;

elseif strcmp(mode_sol, 'efficiency')

    % Compute beamformer target
    [tau, h_fig_tau] = bf_grq_tau_tgt(A, R, freq, tau_tgt, varargin_tau_tgt{:});

    % Compute beamformer weights
    for m = 1:M
        W(:, m) = bf_grq_mecd(A(:,:,m), A(:,:,m), R(:,:,m), tau(m));
    end

elseif strcmp(mode_sol, 'efficiency-nowt')
    
    assert(~isempty(wt), 'Requires non-empty options.wt');

    % Compute beamformer weights for non-weighted GRQ constraints, and weighted efficiency covariance:
    % Maximize frequency weighted W'*A*W s.t. W'*W = 1, GRQ(A_nowt, R_nowt, W) = tau constraints

    % Compute beamformer target
    [tau, h_fig_tau] = bf_grq_tau_tgt(A_nowt, R_nowt, freq, tau_tgt, varargin_tau_tgt{:});

    % Compute beamformer weights
    for m = 1:M
        W(:, m) = bf_grq_mecd(A(:,:,m), A_nowt(:,:,m), R_nowt(:,:,m), tau(m));
    end

    % Revert A, R for display
    A = A_nowt;
    R = R_nowt;

else

    error('Unsupported mode');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_figs = [];
if options.enable_disp && coder.target("MATLAB")
    
    h_figs{1} = h_fig_tau;

    varargin_bf_plt_sh = {'center', options.center, 'wavnum', options.wavnum, 'dB_lim', options.disp_dB_lim, 'title_name', ['GRQ Solver (', mode_sol, ')'], ...
        'disp_freq', options.disp_freq, 'stats_DI_theta_phi', [options.theta_steer, options.phi_steer]};

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
        title_name = ['Solver (', mode_sol, '), modes (A ', mode_A, ', R ', mode_R, '), direction (\theta = ', num2str(rad2deg(options.theta_steer)), ', \phi = ', num2str(rad2deg(options.phi_steer)), ')'];
        h_fig_stats = bf_plt_sh(C_arr, false, W, freq, 'stats', varargin_bf_plt_sh{:}, 'title_name', title_name, ...
            'stats_GRQ_A', A, 'stats_GRQ_R', R, ...
            'DI_dB_lim', options.disp_DI_dB_lim, 'EP_dB_lim', options.disp_EP_dB_lim);
        h_figs{end+1} = h_fig_stats;
    end

end
