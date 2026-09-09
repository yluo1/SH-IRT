function d = sh_cdf_theta_diff_sq(C_pdf, D_pdf, mode, options)
%Compute integrated square difference of marginal cumulative distribution functions 
%of spherical harmonic expansion probability density function:
%CDF of co-latitude theta marginalized over azimuth

%Compute int_{\theta=0}^{\pi} | sh_cdf_theta(C_pdf, \theta) - sh_cdf_theta(D_pdf, \theta) |^2  * sin(theta) d \theta

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:                 [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%D:                 [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%mode:              String, estimation method  {'quadrature', 'analytic', 'analytic_unroll'}

%options:           struct
%options.N_quad:    Number of quadrature points

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

% d_quadrature = sh_cdf_theta_diff_sq(C_pdf_theta_0, C_pdf_theta_half_pi, 'quadrature')
% d_analytic = sh_cdf_theta_diff_sq(C_pdf_theta_0, C_pdf_theta_half_pi, 'analytic')
% norm(d_quadrature - d_analytic)

% d_quadrature = sh_cdf_theta_diff_sq(C_pdf_theta_0, C_pdf_theta_pi, 'quadrature')
% d_analytic = sh_cdf_theta_diff_sq(C_pdf_theta_0, C_pdf_theta_pi, 'analytic')
% norm(d_quadrature - d_analytic)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare across modes for random function

%rng(1561);
%max_odr = 6;
%is_real = false;
%C_pdf = sh_nrm(sh_msq(sh_rand(max_odr, 1, is_real), is_real), 'Sum');
%D_pdf = sh_nrm(sh_msq(sh_rand(max_odr, 1, is_real), is_real), 'Sum');

%d_quadrature = sh_cdf_theta_diff_sq(C_pdf, D_pdf, 'quadrature')
%d_analytic = sh_cdf_theta_diff_sq(C_pdf, D_pdf, 'analytic')
%d_analytic_unroll  = sh_cdf_theta_diff_sq(C_pdf, D_pdf, 'analytic_unroll')

%err_quadrature = norm(d_analytic - d_quadrature)
%err_analytic_unroll = norm(d_analytic - d_analytic_unroll)

arguments
    C_pdf (:,:) double {coder.mustBeComplex} = complex(0);
    D_pdf (:,:) double {coder.mustBeComplex} = complex(0);
    mode (1,:) char {mustBeMember(mode, {'quadrature', 'analytic', 'analytic_unroll'})} = 'analytic';

    options.N_quad (1,1) double {mustBeInteger, mustBePositive} = 1001;
end

[P, M] = size(C_pdf);
P = sqrt(P) - 1;
assert(P == floor(P), 'Invalid size C');

assert(isequal(size(C_pdf), size(D_pdf)), 'C_pdf, D_pdf size mismatch');

if strcmp(mode, 'quadrature')

    N = options.N_quad;
    theta = linspace(0, pi, N)';
    d_theta = pi / (N - 1);

    y = sh_cdf_theta(C_pdf - D_pdf, theta); %[N x M]
    d = d_theta * sin(theta') * (y.^2);

    
elseif strcmp(mode, 'analytic')
    
    E_pdf = C_pdf - D_pdf;

    P_acc = zeros(P+2, M); %Accumulator of coefficients per Legendre polynomial

    for l = 0:P

        idx = (l + 1)^2 - l;

        if l == 0

            coeff_l = sqrt(pi) * real(E_pdf(idx, :));
            P_acc(1, :) = P_acc(1, :) + coeff_l;
            P_acc(2, :) = P_acc(2, :) - coeff_l;

        else %l > 0

            coeff_l = sqrt(pi / (2 * l + 1)) * real(E_pdf(idx, :));
            P_acc(l, :)   = P_acc(l, :)   + coeff_l;
            P_acc(l+2, :) = P_acc(l+2, :) - coeff_l;

        end
    end

    d = (2 ./ ( 2 * (0:(P+1)) + 1 )) * (P_acc.^2);

elseif strcmp(mode, 'analytic_unroll')

    E_pdf = C_pdf - D_pdf;
    
    P_acc = zeros(P+2, M); %Accumulator of coefficients per Legendre polynomial

    for l = 0:(P+1)
        
        if l == 0

            if P == 0
                
                idx_lo = (l + 1)^2 - l;    
                P_acc(1, :) = sqrt(pi) * real(E_pdf(idx_lo, :));                           

            else

                idx_lo = (l + 1)^2 - l;
                idx_hi = ((l + 1) + 1)^2 - (l + 1);
    
                P_acc(1, :) = sqrt(pi) * real(E_pdf(idx_lo, :)) ...
                            + sqrt(pi / (2 + 1)) * real(E_pdf(idx_hi, :));
            end

        elseif l >= P

            idx_lo = ((l - 1) + 1)^2 - (l - 1);
            
            P_acc(l+1, :) = - sqrt(pi / (2 * (l - 1) + 1)) * real(E_pdf(idx_lo, :));                      

        else

            idx_lo = ((l - 1) + 1)^2 - (l - 1);
            idx_hi = ((l + 1) + 1)^2 - (l + 1);

            P_acc(l+1, :) = - sqrt(pi / (2 * (l - 1) + 1)) * real(E_pdf(idx_lo, :)) ...
                            + sqrt(pi / (2 * (l + 1) + 1)) * real(E_pdf(idx_hi, :));         

        end
    end

    d = (2 ./ ( 2 * (0:(P+1)) + 1 )) * (P_acc.^2);

else
    error('Unknown mode');
end
