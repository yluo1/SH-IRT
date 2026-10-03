function plot_sh_fit_msq_example(max_odr_half, num_func, options)
%Plot sh_fit_msq.m modes

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%max_odr_half:      Component max-order expansion in reference field  
%num_func:          Number of sum-of-magnitude square component functions

%options:           struct
%options.rseed:     Randomized seed
%options.dB_lim:    [1 x 2] Display dB range [min, max]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Plot fits for reference non-negative fields of varying
%number of sum-of-magnitude square components 

%plot_sh_fit_msq_example(4, 1);
%plot_sh_fit_msq_example(4, 3);

%plot_sh_fit_msq_example(6, 1);
%plot_sh_fit_msq_example(6, 3);

arguments
    max_odr_half (1,1) double {mustBeNonnegative, mustBeInteger} = 3;
    num_func (1,1) double {mustBePositive, mustBeInteger} = 1;

    options.rseed {mustBeNonnegative, mustBeInteger} = 441;
    options.dB_lim (1,2) double = [0, 40];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate sample sum-of-magnitude square function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(options.rseed);
max_odr = 2 * max_odr_half;
N_C = (max_odr + 1)^2;
N_pts = N_C;

is_real = false;

% Generate reference field
C_refs = sh_rand(max_odr_half, num_func, is_real);
D_ref = sum(sh_msq(C_refs, is_real), 2);

% Sample Fibonnaci uniform points over sphere
[theta, phi] = sh_grd_fib(N_pts);  

% Sample from field, add random noise
X = real(sh_dec(D_ref, theta, phi, is_real));
X = X + randn(N_pts, 1) * 5e-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit non-negative SH expansions to noisy observations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[D_SOMS, C_SOMS]    = sh_fit_msq(X, theta, phi, max_odr, is_real, 'SOMS');
[D_MS, C_MS]        = sh_fit_msq(X, theta, phi, max_odr, is_real, 'MS', 'C0', sh_rand(max_odr_half, 20, is_real) );

B = sh_rand(max_odr_half, N_C, is_real);
[D_MOMS, C_MOMS]    = sh_fit_msq(X, theta, phi, max_odr, is_real, 'MOMS', 'B', B);
[D_MP, C_MP]    = sh_fit_msq(X, theta, phi, max_odr, is_real, 'MP', 'B', B, 'B0', rand(N_C, 20));

err_SOMS    = err_NSHMSQ(D_ref, D_SOMS)
err_MS      = err_NSHMSQ(D_ref, D_MS)
err_MOMS    = err_NSHMSQ(D_ref, D_MOMS)
err_MP      = err_NSHMSQ(D_ref, D_MP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varargin = {'disp_phase', false, 'dB_lim', options.dB_lim };
h_ref   = sh_plt(D_ref, 'mercator', is_real, 'title_name_override', 'Reference', varargin{:}, 'disp_theta_phi', [theta, phi], 'disp_theta_phi_markersize', 8);

h_SOMS  = sh_plt(D_SOMS, 'mercator', is_real, 'title_name_override', ['Sum-of-Magnitude Square Fit NMSQE: ', num2str(err_SOMS)], varargin{:});
h_MS    = sh_plt(D_MS, 'mercator', is_real, 'title_name_override', ['Magnitude Square Fit NMSQE: ', num2str(err_MS)], varargin{:});

h_MOMS  = sh_plt(D_MOMS, 'mercator', is_real, 'title_name_override', ['Mix-of-Magnitude Square Fit NMSQE: ', num2str(err_MOMS)], varargin{:});
h_MP    = sh_plt(D_MP, 'mercator', is_real, 'title_name_override', ['Mixture Power Fit NMSQE: ', num2str(err_MP)], varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_dir = 'figs/figs_sh/fit';
if ~isfolder(out_dir)
    mkdir(out_dir);
end

exportgraphics(h_ref{1}, fullfile(out_dir,  ['msq_ref_', num2str(max_odr_half), '_', num2str(num_func), '.png'] ));
exportgraphics(h_SOMS{1}, fullfile(out_dir, ['msq_soms_', num2str(max_odr_half), '_', num2str(num_func), '.png'] ));
exportgraphics(h_MS{1}, fullfile(out_dir,   ['msq_ms_', num2str(max_odr_half), '_', num2str(num_func), '.png'] ));
exportgraphics(h_MOMS{1}, fullfile(out_dir, ['msq_moms_', num2str(max_odr_half), '_', num2str(num_func), '.png'] ));
exportgraphics(h_MP{1}, fullfile(out_dir,   ['msq_mp_', num2str(max_odr_half), '_', num2str(num_func), '.png'] ));