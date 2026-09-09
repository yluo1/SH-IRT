function d = sh_pdf_dist(C_pdf, D_pdf, mode, is_real, options)
%Compute Wasserstein distance metrics between C_pdf to D_pdf,

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C_pdf:     [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%D_pdf:     [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%mode:      String, compute method {'Wasserstein2'}
%               'SlicedWasserstein2':   Approximate Wasserstein p=2:  W_2^2(C_pdf, D_pdf)  
%                                       via spherical sliced projections onto semi-circle given by 
%                                       int_{\Omega} sh_cdf_theta_diff_sq(R(\Omega) * C_pdf, R(\Omega) * D_pdf) d \Omega
%                                       for SH rotations onto uniform directions \Omega over the sphere 
%               'KLDiv':                Kullback-Lieber Divergence(C_pdf | D_pdf)
%                                       int_{\Omega} f_C(\Omega) * log(f_C(\Omega) / f_D(\Omega)) d\Omega 

%is_real:       Logical, if true, evaluate real SH

%options                   struct
%options.SW_NC_fac:        Oversample number of projections by this non-negative factor
%options.SW_N_CDF_quad:    Number of quadrature points for sampling inverse CDF
%options.SW_mode_quad:     String, quadrature method  {'uniform', 'nonuni_C', 'nonuni_D', 'interp1'}, see sh_cdf_inv_theta_diff_sq.m
%options.SW_mode_slice_wt: String, weighting method per slice {'uniform', 'exp'}
%                               'uniform':      1
%                               'exp':          exp(sh_cdf_theta_diff_sq(R(\Omega) * C_pdf, R(\Omega) * D_pdf))
%                                               Nguyen, K., & Ho, N. (2023). Energy-based sliced wasserstein distance. Advances in Neural Information Processing Systems, 36, 18046-18075.

%options.KLD_NC_fac:        Oversample number of points to sample on sphere by this non-negative factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%d:          [1 x M] distance metric

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare Wasserstein p=2 distance rotation invariance

% max_odr = 3;
% theta = 0;
% phi = 0;
% ell = 0.5;
% is_real = false;
% 
% C_pdf_theta_0 = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, 0, phi, ell, is_real) ), 'Sum');
% sh_plt(C_pdf_theta_0, 'mercator', is_real);
% 
% C_pdf_theta_half_pi = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, pi/2, phi, ell, is_real) ), 'Sum');
% sh_plt(C_pdf_theta_half_pi, 'mercator', is_real);
%
% C_pdf_theta_half_pi_phi_pi = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, pi/2, pi, ell, is_real) ), 'Sum');
% sh_plt(C_pdf_theta_half_pi_phi_pi, 'mercator', is_real);
%
% C_pdf_theta_pi = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, pi, phi, ell, is_real) ), 'Sum');
% sh_plt(C_pdf_theta_pi, 'mercator', is_real);


%d_theta_0_theta_pi = sh_pdf_dist(C_pdf_theta_0, C_pdf_theta_pi, 'SlicedWasserstein2', is_real)
%d_theta_half_pi_theta_half_pi_phi_pi = sh_pdf_dist(C_pdf_theta_half_pi, C_pdf_theta_half_pi_phi_pi, 'SlicedWasserstein2', is_real)
%d_theta_0_theta_half_pi = sh_pdf_dist(C_pdf_theta_0, C_pdf_theta_half_pi, 'SlicedWasserstein2', is_real)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Rotated RBF density kernels

% max_odr = 4;
% is_real = false;
% theta = pi/2;
% phi = 0;
% ell = 0.5;
% M = 50;
% SW_mode_slice_wt = 'uniform';
% %SW_mode_slice_wt = 'exp';
% C_pdf_ref = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, theta, phi,  ell, is_real) ), 'Sum');
% phi_list = linspace(0, 2*pi, M)';
% C_pdf = zeros((2 * max_odr + 1)^2, M);
% for m = 1:M
%    C_pdf(:, m) = sh_rot(C_pdf_ref, [phi_list(m), 0, 0], is_real);
% end
% d_SW = sh_pdf_dist(repmat(C_pdf_ref, [1, M]), C_pdf, 'SlicedWasserstein2', is_real, 'SW_mode_slice_wt', SW_mode_slice_wt);
% d_KLD = sh_pdf_dist(repmat(C_pdf_ref, [1, M]), C_pdf, 'KLDiv', is_real);
% d_KLD_rev = sh_pdf_dist(C_pdf, repmat(C_pdf_ref, [1, M]), 'KLDiv', is_real);
% 
% fontsize = 16;
% figure;
% plot(rad2deg(phi_list), d_SW, 'o-', rad2deg(phi_list), d_KLD, '*-', rad2deg(phi_list), d_KLD_rev, 's-', 'linewidth', 1.5); grid on; axis tight;
% xlabel('Azimuth \phi (Degrees)', 'fontsize', fontsize);
% ylabel('Distance', 'fontsize', fontsize);
% title('Probability Density Function Distances', 'fontsize', fontsize + 1);
% set(gca, 'fontsize', fontsize - 1);
% h_lg = legend('SSW($f_C, R f_C$)', 'KLD($f_C | R f_C$)', 'KLD($R f_C | f_C$)', 'location', 'best', 'interpreter', 'latex');
% set(h_lg, 'fontsize', fontsize - 1);

% sh_plt(C_pdf_ref, 'mercator', is_real);

arguments
    C_pdf (:,:) double {coder.mustBeComplex} = complex(0);
    D_pdf (:,:) double {coder.mustBeComplex} = complex(0);
    mode (1,:) char {mustBeMember(mode, {'SlicedWasserstein2', 'KLDiv'})} = 'SlicedWasserstein2';
    is_real (1,1) logical = false;
    
    options.SW_NC_fac (1,1) double {mustBeNonnegative, mustBeInteger} = 8;
    options.SW_N_CDF_quad (1,1) double {mustBeNonnegative, mustBeInteger} = 251;    
    options.SW_mode_quad (1,:) char {mustBeMember(options.SW_mode_quad,  {'uniform', 'nonuni_C', 'nonuni_D', 'interp1'})} = 'interp1';
    options.SW_mode_slice_wt (1,:) char {mustBeMember(options.SW_mode_slice_wt, {'uniform', 'exp'})} = 'uniform';

    options.KLD_NC_fac (1,1) double {mustBeNonnegative, mustBeInteger} = 32;
end

assert(sh_pdf_check(C_pdf, is_real), 'C_pdf is invalid density');
assert(sh_pdf_check(D_pdf, is_real), 'D_pdf is invalid density');
assert(isequal(size(C_pdf), size(D_pdf)), 'C_pdf and D_pdf size mismatch');

[N_C, M] = size(C_pdf);
P = sqrt(N_C) - 1;
assert(P - floor(P) == 0, 'C_pdf invalid size');

d = zeros([1, M]);

if strcmp(mode, 'SlicedWasserstein2')
    
    N   = options.SW_NC_fac * N_C; %Number of projections

    rot_intrinsic = false;

    [theta, phi] = sh_fib(N);
    psi = zeros(N, 1);

    if ~is_real
        %Convert from complex to real
        C_pdf = sh_cpx2re(C_pdf);
        D_pdf = sh_cpx2re(D_pdf);        
    end

    %Iterate over projections
    total_wt = zeros([1, M]);
    for n = 1:N
        %Compute rotation matrices, real
        C_pdf_rot = sh_rot(C_pdf, -[phi(n), theta(n), psi(n)], true, rot_intrinsic);
        D_pdf_rot = sh_rot(D_pdf, -[phi(n), theta(n), psi(n)], true, rot_intrinsic);

        %d_n = sh_cdf_theta_diff_sq(C_pdf_rot, D_pdf_rot, 'analytic');
        d_n = sh_cdf_inv_theta_diff_sq(C_pdf_rot, D_pdf_rot, options.SW_mode_quad, options.SW_N_CDF_quad);
        
        if strcmp(options.SW_mode_slice_wt, 'uniform')

            total_wt = total_wt + 1/N;

        elseif strcmp(options.SW_mode_slice_wt, 'exp')

            wt_n = exp(d_n);
            d_n = d_n .* wt_n;
            total_wt = total_wt + wt_n;
        
        else
            error('Unsupported options.SW_mode_slice_wt')
        end
        d = d + d_n;
    end
    d = d ./ total_wt;
    d = d / N;

elseif strcmp(mode, 'KLDiv')
    
    N   = options.KLD_NC_fac * N_C; %Number of points on sphere
    [theta, phi] = sh_fib(N);

    f_C = abs(sh_dec(C_pdf, theta, phi));
    f_D = abs(sh_dec(D_pdf, theta, phi));

%    d = mean(f_C .* log(f_C ./ f_D), 1) * (4 * pi);
     d = mean(f_C .* (log(f_C) -  log(f_D)), 1) * (4 * pi);
      
else
    error('Unknown mode');
end

