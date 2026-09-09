function [y, dy] = sh_cdf_phi_cond(C_pdf, theta, phi)
%Compute conditional cumulative distribution function 
%of spherical harmonic expansion probability density function
%of azimuth phi given co-latitude theta

%Compute z(\phi_*)   = \sum_{l=0}^P \sum_{m=-l}^{l} sqrt((2*l + 1)*(l-m)!/(4*pi*(l+m)!)) C_l^m P_l^m(\cos(\theta_*)) \int_{\phi = 0}^{\phi_*} exp(1i * m * phi)
%                    = \sum_{l=0}^P \sum_{m=-l}^{l} sqrt((2*l + 1)*(l-m)!/(4*pi*(l+m)!)) C_l^m P_l^m(\cos(\theta_*)) (sin(m * phi_*))/ m
%For valid PDFs, replace exp(1i * m * phi) with cos(m * phi) as the imaginary components sum to 0

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:                 [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%theta:             [1 x 1] Co-latitude (radians), (0, pi)
%phi:               [N x 1] Azimuth (radians), (0, 2 * pi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%y:                 [N x M] cumulative distribution
%dy:                [N x M] derivative of CDF w.r.t. phi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compute conditional cumulative distribution function for sample SH expansion of 
%squared exponential kernel over various co-latitude and azimuth

% rng(234);
% max_odr = 10;
% ell = 0.5;
% is_real = false;
% 
% C_pdf_theta_0 = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, 0, 0, ell, is_real) ), 'Sum');
% sh_plt(C_pdf_theta_0, 'mercator', is_real);
% 
% C_pdf_theta_half_pi = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, pi/2, 0, ell, is_real) ), 'Sum');
% sh_plt(C_pdf_theta_half_pi, 'mercator', is_real);
% 
% C_pdf_theta_half_pi_rot = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, pi/2, pi, ell, is_real) ), 'Sum');
% sh_plt(C_pdf_theta_half_pi_rot, 'mercator', is_real);
% 
% phi = linspace(0, 2*pi, 100)';
% theta_half_pi = pi/2;
% y_theta_0 = sh_cdf_phi_cond(C_pdf_theta_0, theta_half_pi, phi);
% y_theta_half_pi = sh_cdf_phi_cond(C_pdf_theta_half_pi, theta_half_pi, phi);
% y_theta_half_pi_rot = sh_cdf_phi_cond(C_pdf_theta_half_pi_rot, theta_half_pi, phi);

% fontsize = 14;
% figure; plot(phi, y_theta_0, phi, y_theta_half_pi, phi, y_theta_half_pi_rot, 'linewidth', 1.5); grid on; axis tight;
% xlabel('Azimuth \phi', 'fontsize', fontsize);
% ylabel('CDF(\phi | \theta = \pi/2)', 'fontsize', fontsize);
% title('Spherical Harmonic Expansion Cumulative Distribution', 'fontsize', fontsize +1);
% set(gca, 'fontsize', fontsize - 1);
% h_lg = legend('Exp $\theta = 0$', 'Exp $\theta = \pi/2$, $\phi = 0$', 'Exp $\theta = \pi/2$, $\phi = \pi$', 'location', 'best', 'interpreter', 'latex');
% set(h_lg, 'fontsize', fontsize - 1);

arguments
     C_pdf (:,:) double {coder.mustBeComplex} = complex(0);
     theta (1,1) double {mustBeNonnegative} = [0];
     phi (:,1) double {mustBeNonnegative} = [0];
end

[P, M] = size(C_pdf);
P = sqrt(P) - 1;
assert(P == floor(P), 'Invalid size C');

N = numel(phi);

cos_theta = cos(theta);

%Compute CDF
y = complex(zeros([N, M]));
z = zeros([1, M]); %Normalizing term
for l = 0:P
    leg_l_half = asc_legendre(l, cos_theta);
    if l == 0
        leg_l = leg_l_half;
    else
        idx_m = (1:l)';
        leg_l = [flipud(leg_l_half(2:end) .* (-1).^idx_m .* factorial(l - idx_m) ./ factorial(l + idx_m)); leg_l_half];
    end

    for m = -l:l
        idx = l^2 + l + m + 1;
        if m == 0

          % y = y + sqrt((2*l + 1)*factorial(l - m) / (4*pi*factorial(l + m))) ...
          %      * leg_l(m + l + 1) * phi ...
          %      * (C_pdf(idx, :));
          % 
          % z = z + sqrt((2*l + 1)*factorial(l - m) / (4*pi*factorial(l + m))) ...
          %      * leg_l(m + l + 1) * (2*pi) ...
          %      * (C_pdf(idx, :));

          fac = sqrt((2*l + 1)*factorial(l - m) / (4*pi*factorial(l + m))) ...
               * leg_l(m + l + 1) * real(C_pdf(idx, :));

          y = y + fac * phi;
          z = z + fac * 2 * pi;

        else

           % y = y + sqrt((2*l + 1)*factorial(l - m) / (4*pi*factorial(l + m))) ...
           %     * leg_l(m + l + 1) * (sin(m * phi) - 1i * cos(m * phi) + 1i) / m ...
           %     * (C_pdf(idx, :));
           % 
           % z = z + sqrt((2*l + 1)*factorial(l - m) / (4*pi*factorial(l + m))) ...
           %     * leg_l(m + l + 1) * (sin(m * 2*pi) - 1i * cos(m * 2*pi) + 1i) / m ...
           %     * (C_pdf(idx, :));

            fac = sqrt((2*l + 1)*factorial(l - m) / (4*pi*factorial(l + m))) ...
               * leg_l(m + l + 1) / m * (C_pdf(idx, :));

            y = y + fac * (sin(m * phi) - 1i * cos(m * phi) + 1i);

        end
    end    
end
y = real(bsxfun(@rdivide, y, z)); %Normalize

if nargout > 1  % Compute derivative of CDF
    dy = real(sh_dec(C_pdf, theta * ones(N, 1), phi, false));
    dy = real(bsxfun(@rdivide, dy, z)); %Normalize
end