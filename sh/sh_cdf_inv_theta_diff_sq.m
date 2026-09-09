function d = sh_cdf_inv_theta_diff_sq(C_pdf, D_pdf, mode, N)
%Compute integrated square difference of marginal cumulative inverse distribution functions 
%of spherical harmonic expansion probability density function:
%CDF of co-latitude theta marginalized over azimuth

%Compute int_{u=0}^{1} | sh_cdf_inv_theta(C_pdf, u) - sh_cdf_inv_theta(D_pdf, u) |^2 du

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:                 [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%D:                 [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)

%mode:              String, quadrature method  {'uniform', 'nonuni_C', 'nonuni_D', 'interp1'}
%                       'uniform':          Quadrature by uniform spaced u sampling of sh_cdf_inv_theta (Reference)
%                       'nonuni_C':         Quadrature by non-uniform spaced u = sh_cdf_theta(C_pdf, theta) for uniform theta
%                       'nonuni_D':         Quadrature by non-uniform spaced u = sh_cdf_theta(D_pdf, theta) for uniform theta
%                       'interp1':          Quadrature by uniform spaced u linear interpolation of 
%                                           piece-wise linear approximation of sh_cdf_inv_theta (Fastest)

%N:                 Number of quadrature points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%d:                 [1 x M] Integrated square difference

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:      CDF distance between RBF kernels

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
% C_pdf_theta_pi = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, pi, phi, ell, is_real) ), 'Sum');
% sh_plt(C_pdf_theta_pi, 'mercator', is_real);
% 
% theta_s = linspace(0, pi, 100)';
% y_theta_0 = sh_cdf_theta(C_pdf_theta_0, theta_s);
% y_theta_half_pi = sh_cdf_theta(C_pdf_theta_half_pi, theta_s);
% y_theta_pi = sh_cdf_theta(C_pdf_theta_pi, theta_s);
%
% fontsize = 14;
% figure; plot(theta_s, y_theta_0, theta_s, y_theta_half_pi, theta_s, y_theta_pi, 'linewidth', 1.5); grid on; axis tight; 
% xlabel('Co-latitude \theta', 'fontsize', fontsize);
% ylabel('CDF(\theta)', 'fontsize', fontsize);
% title('Spherical Harmonic Expansion Cumulative Distribution', 'fontsize', fontsize +1);
% set(gca, 'fontsize', fontsize - 1);
% h_lg = legend('Exp $\theta = 0$', 'Exp $\theta = \pi/2$', 'Exp $\theta = \pi$', 'location', 'best', 'interpreter', 'latex');
% set(h_lg, 'fontsize', fontsize - 1);

% d_uniform = sh_cdf_inv_theta_diff_sq(C_pdf_theta_0, C_pdf_theta_half_pi, 'uniform')
% d_uniform = sh_cdf_inv_theta_diff_sq(C_pdf_theta_0, C_pdf_theta_pi, 'uniform')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare across modes for random function

%rng(1561);
%max_odr = 6;
%is_real = false;
%C_pdf = sh_nrm(sh_msq(sh_rand(max_odr, 1, is_real), is_real), 'Sum');
%D_pdf = sh_nrm(sh_msq(sh_rand(max_odr, 1, is_real), is_real), 'Sum');

%tic; d_uniform = sh_cdf_inv_theta_diff_sq(C_pdf, D_pdf, 'uniform')
%toc_uniform = toc;
%tic; d_nonuni_C = sh_cdf_inv_theta_diff_sq(C_pdf, D_pdf, 'nonuni_C')
%toc_nonuni_C = toc;
%tic; d_nonuni_D = sh_cdf_inv_theta_diff_sq(C_pdf, D_pdf, 'nonuni_D')
%toc_nonuni_D = toc;
%tic; d_interp1 = sh_cdf_inv_theta_diff_sq(C_pdf, D_pdf, 'interp1')
%toc_interp1 = toc;

%toc_uniform
%toc_nonuni_C
%toc_nonuni_D
%toc_interp1
%err_nonuni_C = norm(d_uniform - d_nonuni_C)
%err_nonuni_D = norm(d_uniform - d_nonuni_D)
%err_interp1 = norm(d_uniform - d_interp1)

arguments
    C_pdf (:,:) double {coder.mustBeComplex} = complex(0);
    D_pdf (:,:) double {coder.mustBeComplex} = complex(0);
    mode (1,:) char {mustBeMember(mode, {'uniform', 'nonuni_C', 'nonuni_D', 'interp1'})} = 'interp1';
    N (1,1) double {mustBeInteger, mustBePositive} = 251;
end

[P, M] = size(C_pdf);
P = sqrt(P) - 1;
assert(P == floor(P), 'Invalid size C');

assert(isequal(size(C_pdf), size(D_pdf)), 'C_pdf, D_pdf size mismatch');

if strcmp(mode, 'uniform') %uniform spaced u

    u = linspace(0, 1, N)';
    d_u = 1 / (N - 1);

    theta_C = sh_cdf_inv_theta(C_pdf, u); %[N x M]
    theta_D = sh_cdf_inv_theta(D_pdf, u); %[N x M]

    d = d_u * ones(1, N) * ((theta_C - theta_D).^2);
    
elseif strcmp(mode, 'nonuni_C') %non-uniform spaced u = sh_cdf_theta(C_pdf, theta) for uniform spaced theta

    theta = linspace(0, pi, N)';

    u_C = sh_cdf_theta(C_pdf, theta); %[N x M]
    d = zeros(1, M);
    for m = 1:M
        theta_D = sh_cdf_inv_theta(D_pdf, u_C); 
        d(m) = diff(u_C(:, m))' * ((theta(2:end) - theta_D(2:end)).^2);
    end

elseif strcmp(mode, 'nonuni_D') %non-uniform spaced u = sh_cdf_theta(D_pdf, theta) for uniform spaced theta

    theta = linspace(0, pi, N)';

    u_D = sh_cdf_theta(D_pdf, theta); %[N x M]
    d = zeros(1, M);
    for m = 1:M
        theta_C = sh_cdf_inv_theta(C_pdf, u_D); 
        d(m) = diff(u_D(:, m))' * ((theta(2:end) - theta_C(2:end)).^2);
    end

elseif strcmp(mode, 'interp1') %Approximate sh_cdf_inv_theta via linear-interpolations 
                               %of u_n = sh_cdf_theta(*, theta_n)
    
    theta = linspace(0, pi, N)';
    u_C = sh_cdf_theta(C_pdf, theta); %[N x M]
    u_D = sh_cdf_theta(D_pdf, theta); %[N x M]

    %Replace end-points
    u_C(1, :) = 0;
    u_C(end, :) = 1;
    u_D(1, :) = 0;
    u_D(end, :) = 1;

    u = linspace(0, 1, N)';
    d_u = 1 / (N - 1);

    d = zeros(1, M);
    for m = 1:M
        d(m) = d_u * ones(1, N) * (interp1(u_C(:, m), theta, u, 'linear') - interp1(u_D(:, m), theta, u, 'linear')).^2;
    end
    
else
    error('Unknown mode');
end
