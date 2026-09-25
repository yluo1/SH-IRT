function tst_sh_fit_msq(max_odr_half, is_real, is_pdf, N, options)
%Test sh_fit_msq.m

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%max_odr_half:  max-order for SH C (half of D)
%is_real:       Logical, if true, evaluate real SH     
%is_pdf:        Logical, if true, constraint fits to be density functions
%N:             Number of sample points

%options:               struct
%options.rseed:         Random seed
%options.enable_disp:   Logical, if true, plot target and fitted expansions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:

%tst_sh_fit_msq(4, true, true, 100);
%tst_sh_fit_msq(4, true, false, 100);
%tst_sh_fit_msq(4, false, true, 100);
%tst_sh_fit_msq(4, false, false, 100);

arguments
    max_odr_half (1,1) double {mustBeNonnegative, mustBeInteger} = 4;
    is_real (1,1) logical = true;
    is_pdf (1,1) logical = true;
    N (1,1) double {mustBePositive} = 100;
    options.rseed (1,1) double {mustBeNonnegative, mustBeInteger} = 23541;
    
    options.enable_disp (1,1) logical = true
end

rng(options.rseed);

max_odr = 2 * max_odr_half;

C_ref = sh_rand(max_odr_half, 1, is_real);
if is_pdf
    C_ref = C_ref / sqrt(C_ref'*C_ref);
end
D_ref = sh_msq(C_ref, is_real);

[theta, phi] = sh_fib(N);
X = real(sh_dec(D_ref, theta, phi, is_real));

K = (max_odr_half + 1)^2;
B = sh_rand(max_odr_half, K, is_real);
C0 = sh_rand(max_odr_half, K, is_real);
if is_real
    B0 = randn(K, K);
else
%    C0 = sh_rand(max_odr_half, K, is_real);
    B0 = randn(K, K) + randn(K, K) * 1i;
end

[D_MS, C_MS, ~]                     = sh_fit_msq(X, theta, phi, max_odr, is_real, 'MS', 'C0', [C0], 'is_pdf', is_pdf);
[D_SOMS, C_SOMS, ~]                 = sh_fit_msq(X, theta, phi, max_odr, is_real, 'SOMS', 'is_pdf', is_pdf);
[D_MOMS, C_MOMS, ~, output_MOMS]    = sh_fit_msq(X, theta, phi, max_odr, is_real, 'MOMS', 'B', [B], 'is_pdf', is_pdf);
[D_MP, C_MP, ~, output_MP]          = sh_fit_msq(X, theta, phi, max_odr, is_real, 'MP', 'B', [B], 'B0', B0, 'is_pdf', is_pdf);

err_MS = err_SHMSQ(D_ref, D_MS)
err_SOMS = err_SHMSQ(D_ref, D_SOMS)
err_MOMS = err_SHMSQ(D_ref, D_MOMS)
err_MP = err_SHMSQ(D_ref, D_MP)

if is_pdf
    is_pdf_MS = sh_pdf_check(D_MS, is_real)
    is_pdf_SOMS = sh_pdf_check(D_SOMS, is_real)
    is_pdf_MOMS = sh_pdf_check(D_MOMS, is_real)
    is_pdf_MP = sh_pdf_check(D_MP, is_real)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if options.enable_disp
    varargin = {'disp_theta_phi', [theta, phi], 'disp_theta_phi_markersize', 24, 'dB_lim', [-80, 20]};

    % Plot C
    % sh_plt(C_ref, 'mercator', is_real, 'title_name', 'Ref. C', varargin{:});
    % sh_plt(C_MS, 'mercator', is_real, 'title_name', 'Magnitude Square C', varargin{:});
    % sh_plt(C_MP, 'mercator', is_real, 'title_name', 'Mixture Power C', varargin{:});
    % sh_plt(C_SOMS, 'mercator', is_real, 'title_name', 'Sum-of-Magnitude Square C', varargin{:});
    % sh_plt(C_MOMS, 'mercator', is_real, 'title_name', 'Mix-of-Magnitude Square C', varargin{:});
    
    % Plot D
    sh_plt(D_ref, 'mercator', is_real, 'title_name', 'Ref. D', varargin{:});
    sh_plt(D_MS, 'mercator', is_real, 'title_name', 'Magnitude Square D', varargin{:});
    sh_plt(D_SOMS, 'mercator', is_real, 'title_name', 'Sum-of-Magnitude Square  D', varargin{:});
    sh_plt(D_MOMS, 'mercator', is_real, 'title_name', 'Mix-of-Magnitude Square D', varargin{:});
    sh_plt(D_MP, 'mercator', is_real, 'title_name', 'Mixture Power D', varargin{:});
end
;