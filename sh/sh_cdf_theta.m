function [y, dy] = sh_cdf_theta(C_pdf, theta)
%Compute marginal cumulative distribution function 
%of spherical harmonic expansion probability density function:
%CDF of co-latitude theta marginalized over azimuth

%Compute y(\theta_*) = \sum_{l=0}^P \sum_{m=-l}^l sqrt((2*l + 1)*(l-m)!/(4*pi*(l+m)!)) C_l^m
%                    * \int_{\theta=0}^{\theta_*} \int_{\phi = 0}^{2 pi}   sin(\theta) P_l^m(\cos(\theta)) e^(1i m phi) d\theta d\phi

%                    = \sum_{l=0}^P sqrt((2*l + 1)/(4*pi)) C_l^0 * ( (P_{l+1}(\cos(\theta_*)) - P_{l-1}(\cos(\theta_*))) -  (P_{l+1}(\cos(0)) - P_{l-1}(\cos(0))) ) / (2l + 1)
%Restrict m = 0 as the marginal follows int_{phi=0}^{2 * pi} e^(1i * m * phi) = 0 for m ~= 0

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:                 [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%theta:             [N x 1] Co-latitude (radians), (0, pi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%y:                 [N x M] cumulative distribution
%dy:                [N x M] derivative of CDF w.r.t. theta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compute cumulative distribution function for sample SH expansion of 
%squared exponential kernel over various co-latitude

% rng(234);
% max_odr = 10;
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

arguments
     C_pdf (:,:) double {coder.mustBeComplex} = complex(0);
     theta (:,1) double {mustBeNonnegative} = [0];
end

[P, M] = size(C_pdf);
P = sqrt(P) - 1;
assert(P == floor(P), 'Invalid size C');

assert(all(theta >= 0) && all(theta <= pi), 'theta must be between [0, pi]');

N = numel(theta);

cos_theta = cos(theta);

%Compute CDF
y = zeros([N, M]);
for l = 0:P
    idx = (l + 1)^2 - l;
    if l == 0
        % y = y + sqrt((2*l + 1) / (4*pi)) ...
        %     * ( (1 - legendre_poly(l+1, cos_theta) ) - (1 - legendre_poly(l+1, ones(N, 1)) )  ) / (2*l + 1) ...
        %     * real(C_pdf(idx, :));  %Requires normalization

        % y = y + sqrt((2*l + 1) / (4*pi) ) ...
        %     * ( 1 - legendre_poly(l+1, cos_theta)  ) / (2*l + 1) ...
        %     * real(C_pdf(idx, :)); %Requires normalization

        y = y + sqrt(pi / (2*l + 1)) ...
            * ( 1 - legendre_poly(l+1, cos_theta)  ) ...
            * real(C_pdf(idx, :)); %No normalization

    else
        % y = y + sqrt((2*l + 1) / (4*pi)) ...
        %     * ( (legendre_poly(l-1, cos_theta) - legendre_poly(l+1, cos_theta) ) - (legendre_poly(l-1, ones(N, 1)) - legendre_poly(l+1, ones(N, 1)) )  ) / (2*l + 1) ...
        %     * real(C_pdf(idx, :));  %Requires normalization

        % y = y + sqrt((2*l + 1) / (4*pi)) ...
        %     * ( (legendre_poly(l-1, cos_theta) - legendre_poly(l+1, cos_theta) ) ) / (2*l + 1) ...
        %     * real(C_pdf(idx, :)); %Requires normalization

        y = y + sqrt(pi / (2*l + 1)) ...
            * ( (legendre_poly(l-1, cos_theta) - legendre_poly(l+1, cos_theta) ) ) ...
            * real(C_pdf(idx, :)); %No normalization
    end
end
%y = y * 2 * pi; %Normalize

if nargout > 1 % Compute derivative of CDF

    dy = zeros([N, M]);
    sin_theta = sin(theta); %[N x 1]
    for l = 0:P
        idx = (l + 1)^2 - l;
%        dy = dy + sqrt((2*l + 1) / (4*pi)) * legendre_poly(l, cos_theta) *  real(C_pdf(idx, :)); %Requires normalization
        dy = dy + sqrt((2*l + 1) * pi) * legendre_poly(l, cos_theta) *  real(C_pdf(idx, :)); %No normalization

    end
    dy = bsxfun(@times, dy, sin_theta);
    %y = dy * 2 * pi; %Normalize 
end
