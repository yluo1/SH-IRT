function [tau, h_fig] = bf_grq_tau_tgt(A, R, freq, tau_tgt, options)
% Generate bounded and smooth generalized Rayleigh quotient targets tau over frequency
% GRQ_lim(1, m) <= tau(m) <= GRQ_lim(2, m) for frequency bins m
% from array independent tau_tgt

% Minimizes squared distance ||tau - tau_tgt||^2  + lambda * ||tau||_{sq accel}^2 (squared acceleration penalty)
% s.t. GRQ_lim(1, m) <= tau(m) <= GRQ_lim(2, m)

% Special cases:
% Input tau_tgt == inf      yields  tau(m) = GRQ_lim(2, m)
% Input tau_tgt == 0        yields  tau(m) = GRQ_lim(1, m)

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%A:         [N x N x M]     Accept covariance matrix
%R:         [N x N x M]     Reject covariance matrix

%freq:      [1 x M]  Frequency in Hz

%tau_tgt:   Scalar          Constant target                           
%           [1 x M]         Target per frequency

%options:                   Struct
%options.U:                 Empty or [N x N_U]  Null space vectors
%                           Disable constraints if empty

%options.lambda:            Scalar, regularization term

%options.enable_disp:       Logical, if true, display targets
%options.disp_dB_lim:       [1 x 2] Display [min, max] dB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%tau:       [1 x M] GRQ targets
%h_fig:     Handle to figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate GRQ targets for phased-array

% Fs = 16000;
% M = 256;
% N = 4;
% radius = 0.1;
% delay_dist = 1;
% mode = 'uni_circ';
% freq = logspace(log10(20), log10(Fs), M + 1); freq = freq(1:end-1);
% wavnum = 2 * pi * freq / 343;
% is_real = false;
% min_eig = 1e-8;

% [C_arr, center] = sh_enc_arr(0, mode, N, radius, delay_dist, freq, 'noncentered', true, 'enable_disp', true, 'disp_dB_lim', [-3, 3]);

% R = sh_cov(C_arr, sh_pdf_preset('FwdCenter', 'Uniform', 0), is_real, 'SHT_center', center, 'SHT_wavnum', wavnum, ...
%     'SHT_oversample', 1, 'SHT_maxsample', 4096, 'min_eig', min_eig);
% A = sh_cov(C_arr, sh_pdf_preset('FwdCenter', 'SqExp', 8, is_real, 'ell', 0.5), is_real, 'SHT_center', center, 'SHT_wavnum', wavnum, ...
%     'SHT_oversample', 1, 'SHT_maxsample', 4096, 'min_eig', min_eig);

%tau_tgt = db2pow(2);
%[tau, h_fig] = bf_grq_tau_tgt(A, R, freq, tau_tgt, 'lambda', 1e4, 'enable_disp', true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate GRQ targets for piston on spherical baffle array

% Fs = 24000;
% M = 256;
% N = 4;
% radius = 0.0254 * 3; % 3 inch
% delay_dist = 1;
% max_odr = 30;
% mode = 'spis_circ';
% freq = logspace(log10(20), log10(Fs), M + 1); freq = freq(1:end-1);
% %freq = linspace(0, Fs, M + 1); freq = freq(1:end-1);
% spis_arc_len = 0.0254;
% is_real = false;
% min_eig = 1e-8;

% C_arr = sh_enc_arr(max_odr, mode, N, radius, delay_dist, freq, 'spis_arc_len', spis_arc_len, 'enable_disp', true, 'disp_dB_lim', [-36, 12]);

% R = sh_cov(C_arr, sh_pdf_preset('FwdCenter', 'Uniform', 0), is_real, ...
%     'SHT_oversample', 1, 'SHT_maxsample', 4096, 'min_eig', min_eig);
% A = sh_cov(C_arr, sh_pdf_preset('FwdCenter', 'SqExp', 8, is_real, 'ell', 0.5), is_real, ...
%     'SHT_oversample', 1, 'SHT_maxsample', 4096, 'min_eig', min_eig);

%tau_tgt = db2pow(5);
%[tau, h_fig] = bf_grq_tau_tgt(A, R, freq, tau_tgt, 'lambda', 1e4, 'enable_disp', true);

arguments
    A (:,:,:) double {coder.mustBeComplex} = complex(1);
    R (:,:,:) double {coder.mustBeComplex} = complex(1);
    freq (1,:) double {mustBeNonnegative} = 0;
    tau_tgt (1,:) {mustBeNonnegative} = db2pow(6);

    options.U (:,:) double {coder.mustBeComplex} = [];

    options.lambda (1,1) {mustBeNonnegative} = 1e4;

    options.enable_disp (1,1) logical = false;
    options.disp_dB_lim (1,2) = [-inf inf];
end

[~, N, M] = size(A);
assert(size(A, 1) == N, 'A must be square');
assert(isequal(size(A), size(R)), 'R size mismatch A');
assert(numel(freq) == M, 'freq mismatch size M of A');

% Compute GRQ limits over frequency
tau_lim = zeros(2, M);
for m = 1:M
    [~, tau_m] = bf_grq_ex( A(:,:,m), R(:,:,m) );
    tau_lim(:, m) = max([tau_m(1); tau_m(2)], 0);  % Truncate
end
tau_lim_dB = pow2db(tau_lim);

tau_tgt_orig = tau_tgt;
if isscalar(tau_tgt) % Constant target
    tau_tgt = ones(1,M) * tau_tgt; 
end
assert(numel(tau_tgt) == M, 'tau_tgt must be scalar or size [1 x M]');
tau_tgt_dB = pow2db(tau_tgt);

% Compute GRQ targets
if isequal(tau_tgt_orig, inf) % Max 

    tau = tau_lim(2,:);

elseif isequal(tau_tgt_orig, 0) % Min 

    tau = tau_lim(1,:);

else

    lambda = options.lambda;

    assert(M >= 5, 'Requires M >= 5');
    H_regu = toeplitz([6, -4, 1, zeros(1, M - 3)]); % Squared acceleration kernel 
    % (assume uniform frequency or log-frequency spacing)
    
    H = eye(M) + lambda * H_regu;
    f = -tau_tgt_dB(:);
 
    %Dirchlet boundary conditions at first two and last two frequencies where
    %the finite-difference kernel intersects past the boundary   
    Aeq = [...
           1, zeros(1, M-1); ...
           0, 1, zeros(1, M-2); ...
           zeros(1, M-2), 1, 0; ...
           zeros(1, M-1), 1; ...
           ];

    beq = [...
           min(max(tau_tgt_dB(1), tau_lim_dB(1,1)), tau_lim_dB(2,1)); ...
           min(max(tau_tgt_dB(2), tau_lim_dB(1,2)), tau_lim_dB(2,2)); ...
           min(max(tau_tgt_dB(M-1), tau_lim_dB(1,M-1)), tau_lim_dB(2,M-1)); ...
           min(max(tau_tgt_dB(M), tau_lim_dB(1,M)), tau_lim_dB(2,M)); ...
           ];

    lb = tau_lim_dB(1, :)';
    ub = tau_lim_dB(2, :)';
    tau_dB = quadprog(H, f, [], [], Aeq, beq, lb, ub); 
    tau = db2pow(tau_dB);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_fig = [];
if options.enable_disp && coder.target("MATLAB")
   
    tau_dB = pow2db(tau);

    fontsize = 16;
    h_fig = figure; 

    semilogx(freq, tau_lim_dB(1,:), 'g--', freq, tau_lim_dB(2,:), 'g-.', ...
       freq, tau_tgt_dB, 'r-',  freq, tau_dB, 'b-.',  'linewidth', 2);
    xlabel('Frequency (Hz)', 'fontsize', fontsize);
    ylabel('Generalized Rayleigh Quotient (dB)', 'fontsize', fontsize);
    grid on; axis tight;
    set(gca, 'fontsize', fontsize - 1);
    h_lg = legend('Upper Bound', 'Lower Bound', 'Input Target', 'Output Target', 'location', 'best');
    set(h_lg, 'fontsize', fontsize - 1);

    title('Generalized Rayleigh Quotient Targets', 'fontsize', fontsize + 1);
    ylim(options.disp_dB_lim);

end
