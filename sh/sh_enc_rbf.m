function C = sh_enc_rbf(mode, max_odr, theta, phi, ell, is_real, enable_disp)
%Spherical harmonic expansion of radial basis functions (RBFs) at spherical coordinates

%RBF function of chordal distances are positive semi-definite and valid kernel functions

%Squared exponential:               exp(-d^2 / (2 * ell^2))
%Matern nu = 5/2:                   (1 + sqrt(5) * d / ell + 5/3 * (d / ell).^2 ) .* exp(-sqrt(5) * d/ell);    
%Matern nu = 3/2:                   (1 + sqrt(3) * d / ell) .* exp(-sqrt(3) * d / ell );    
%Exponential:                       exp(-d / ell)
%Unnormalized sinc:                 sin(ell * d) / (ell * d) 

%d is chordal distance between two vectors u, v on unit sphere:
%d = 2 * sin(|x(u, v)| / 2)  = norm(u - v)

%For unit sphere, x is a angle between u, v:
%x(u, v) = acos(u'v), where v = sph2cart(phi, pi/2 - theta) is fixed, u is free

%Reference:
%Y. Luo, "Spherical harmonic covariance and magnitude function encodings for beamformer design," 
%EURASIP Journal on Audio, Speech, and Music Processing. 2021. 10.1186/s13636-021-00230-7. 

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%mode:          String, RBF function

%               'SqExp':                Squared exponential
%               'SqExpNorm':            Normalized squared exponential
%               'SqExpNonneg':          Squared exponential using squared magnitude expansion
%               'SqExpNonnegNorm':      Normalized squared exponential using squared magnitude expansion
%               'Mat52':                Matern (\nu = 5/2)
%               'Mat32':                Matern (\nu = 3/2)
%               'Exp':                  Exponential
%               'ExpNorm':              Normalized exponential
%               'ExpNonneg':            Exponential using squared magnitude expansion
%               'ExpNonnegNorm':        Normalized exponential using squared magnitude expansion
%               'Sinc':                 Sinc (unnormalized)
%               'SincNum':              Sinc with numerical integration

%max_odr:       Max SH order 
%theta:         [N x 1]  Co-latitude (0, pi)
%phi:           [N x 1]  Azimuth (0, 2 * pi)
%ell:           Scalar, bandwidth parameter, positive

%is_real:       Logical, if true, evaluate real SH

%enable_disp:   Logical, if true, plot function expansion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:             [(max_odr + 1)^2 x N] SH coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_enc_rbf', '-o', 'sh/sh_enc_rbf_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate SH RBF expansions and compare with reference functions

% max_odr = 10;
% theta = pi/2;
% phi = 0;
% ell = 1;
% is_real = false;
% enable_disp = true;
% 
% C_SqExp = sh_enc_rbf('SqExp', max_odr, theta, phi, ell, is_real, enable_disp);
% C_SqExpNonneg = sh_enc_rbf('SqExpNonneg', max_odr, theta, phi, ell, is_real, enable_disp);
% 
% C_Mat52 = sh_enc_rbf('Mat52', max_odr, theta, phi, ell, is_real, enable_disp);
% C_Mat32 = sh_enc_rbf('Mat32', max_odr, theta, phi, ell, is_real, enable_disp);
% 
% C_Exp = sh_enc_rbf('Exp', max_odr, theta, phi, ell, is_real, enable_disp);
% C_ExpNonneg = sh_enc_rbf('ExpNonneg', max_odr, theta, phi, ell, is_real, enable_disp);
% 
% C_Sinc = sh_enc_rbf('Sinc', max_odr, theta, phi, ell / 2, is_real, enable_disp);
% C_SincNum = sh_enc_rbf('SincNum', max_odr, theta, phi, ell / 2, is_real, enable_disp);
% 
% rng(543);
% N_eval = 1000;
% theta_eval = rand(N_eval, 1) * pi; 
% phi_eval   = rand(N_eval, 1) * 2 * pi;
% 
% err_SqExp = norm(sh_dec(C_SqExp, theta_eval, phi_eval, is_real) - srbf_val('SqExp', theta, phi, ell, theta_eval, phi_eval))
% err_Mat52 = norm(sh_dec(C_Mat52, theta_eval, phi_eval, is_real) - srbf_val('Mat52', theta, phi, ell, theta_eval, phi_eval))
% err_Mat32 = norm(sh_dec(C_Mat32, theta_eval, phi_eval, is_real) - srbf_val('Mat32', theta, phi, ell, theta_eval, phi_eval))
% err_Exp   = norm(sh_dec(C_Exp, theta_eval, phi_eval, is_real)   - srbf_val('Exp',   theta, phi, ell, theta_eval, phi_eval))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare Matlab versus Mex

%max_odr = 5;
%theta = pi/2;
%phi = 0;
%ell = 1;
%is_real = false;
%enable_disp = false;

%tic; C_SqExp_mat = sh_enc_rbf('SqExp', max_odr, theta, phi, ell, is_real, enable_disp); toc_SqExp_mat = toc
%tic; C_SqExp_mex = sh_enc_rbf_mex('SqExp', max_odr, theta, phi, ell, is_real, enable_disp); toc_SqExp_mex = toc
%err_SqExp = norm(C_SqExp_mat - C_SqExp_mex)

%tic; C_Mat52_mat = sh_enc_rbf('Mat52', max_odr, theta, phi, ell, is_real, enable_disp); toc_Mat52_mat = toc
%tic; C_Mat52_mex = sh_enc_rbf_mex('Mat52', max_odr, theta, phi, ell, is_real, enable_disp); toc_Mat52_mex = toc
%err_Mat52 = norm(C_Mat52_mat - C_Mat52_mex)

%tic; C_Mat32_mat = sh_enc_rbf('Mat32', max_odr, theta, phi, ell, is_real, enable_disp); toc_Mat32_mat = toc
%tic; C_Mat32_mex = sh_enc_rbf_mex('Mat32', max_odr, theta, phi, ell, is_real, enable_disp); toc_Mat32_mex = toc
%err_Mat32 = norm(C_Mat32_mat - C_Mat32_mex)

%tic; C_Exp_mat = sh_enc_rbf('Exp', max_odr, theta, phi, ell, is_real, enable_disp); toc_Exp_mat = toc
%tic; C_Exp_mex = sh_enc_rbf_mex('Exp', max_odr, theta, phi, ell, is_real, enable_disp); toc_Exp_mex = toc
%err_Exp = norm(C_Exp_mat - C_Exp_mex)

%tic; C_Sinc_mat = sh_enc_rbf('Sinc', max_odr, theta, phi, ell, is_real, enable_disp); toc_Sinc_mat = toc
%tic; C_Sinc_mex = sh_enc_rbf_mex('Sinc', max_odr, theta, phi, ell, is_real, enable_disp); toc_Sinc_mex = toc
%err_Sinc = norm(C_Sinc_mat - C_Sinc_mex)

%tic; C_SincNum_mat = sh_enc_rbf('SincNum', max_odr, theta, phi, ell, is_real, enable_disp); toc_SincNum_mat = toc
%tic; C_SincNum_mex = sh_enc_rbf_mex('SincNum', max_odr, theta, phi, ell, is_real, enable_disp); toc_SincNum_mex = toc
%err_SincNum = norm(C_SincNum_mat - C_SincNum_mex)

arguments
    mode (1,:) char {mustBeMember(mode, {'SqExp', 'SqExpNorm', 'SqExpNonneg', 'SqExpNonnegNorm', 'Mat52', 'Mat32', 'Exp', 'ExpNorm', 'ExpNonneg', 'ExpNonnegNorm', 'Sinc', 'SincNum'})}  = 'SqExp';

    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;

    theta (:,1) double = [0];
    phi   (:,1) double = [0];

    ell (1,1) double {mustBePositive} = 1;

    is_real (1,1) logical = false;

    enable_disp (1,1) logical = false;
    
end

assert(numel(theta) == numel(phi), 'size theta, phi mismatch');
N = numel(theta);

%Compute Legendre polynomial weights as RBF function are zonal w.r.t. z-axis
b = complex(zeros(max_odr + 1, 1));
K = 20;     %Internal number of truncation terms

for l = 0:max_odr

    if strcmp(mode, 'Sinc') %Special case

         b(l+1) = (2 * l + 1)  * pi / (2 * ell) * besselj(l + 1 / 2, ell)^2 ; 

    else %Compute terms

        for k = 0:K
            s = compute_Dnk(l, k) * compute_Ank(l, k, ell, mode);
            b(l+1) = b(l+1) + s;
        end

        b(l+1) = (2 * l + 1) / 2 * b(l+1);
    end
end

%Addition theorem splits Legendre polynomial into SH evaluated at fixed and free angles
%Steer the function by choosing fixed angle to be theta, phi
C = sh_val(max_odr, theta, phi, is_real)';         %Conjugate
for L = 0:max_odr
    idx = (L^2 + 1) : (L+1)^2;
    C(idx, :) = C(idx, :) * b(L+1) * 4 * pi / (2 * L + 1);
end

if contains(mode, 'Nonneg')         %Enforce non-negative function via SH multiplication of the half max order expansion with its conjugate

    idx_half = 1:(floor(max_odr / 2) + 1)^2;

    C = sh_mul( sh_conj(C(idx_half, :)), C(idx_half, :) );
    C = [C; zeros( (max_odr + 1)^2 - size(C, 1), size(C, 2) ) ];

    if contains(mode, 'Norm')       %Normalize
        C = bsxfun(@rdivide, C, sh_int(C, 'Sum'));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if enable_disp && coder.target('MATLAB') 

    x = linspace(0, pi, 512);
    d = 2 * sin(abs(x)/2);

    y_ref = rbf_val(mode, d, ell);
 
    theta_x = theta + x;
    phi_x   = phi * ones(size(x));

    y = real(sh_dec(C, theta_x(:), phi_x(:), is_real))';

    %Compare with reference function
    fontsize = 14;
    figure;
    plot(x, y_ref, 'r-', x, y, 'b-', 'linewidth', 1.5);
    xlabel('Distance d', 'fontsize', fontsize);    
    ylabel('RBF f(d)', 'fontsize', fontsize);
    title(mode, 'fontsize', fontsize + 1);
    grid on; axis tight;
    set(gca, 'fontsize', fontsize - 1);
    h_lg = legend('Reference', 'Expansion', 'location', 'best');
    set(h_lg, 'fontsize', fontsize - 1);

    %Mercator plot
    sh_plt(C, 'mercator', is_real);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dnk] = compute_Dnk(n, k)

nr = 1:n;
kr = 1:k;
Dnk = 4 / pi * exp(sum(log([2*nr, 2*kr - 1, n + kr])) - sum(log([2*nr + 1, kr, 2*(n + kr) + 1])));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ank] = compute_Ank(n, k, p, mode)

if contains(mode, 'SqExp')

    if strcmp(mode, 'SqExpNorm')
        
        Ank = 1/2 * (ank_SqExpNorm(n + 2 * k,  p) - ank_SqExpNorm(n + 2 * k + 2, p) );

    elseif contains(mode, 'SqExpNonneg') 
        
        Ank = 1/2 * (ank_SqExpNonneg(n + 2 * k,  p) - ank_SqExpNonneg(n + 2 * k + 2, p) );

    else

        Ank = 1/2 * (ank_SqExp(n + 2 * k,  p) - ank_SqExp(n + 2 * k + 2, p) );

    end

elseif strcmp(mode, 'Mat52')

    Ank = 1/2 * (ank_Mat52(n + 2 * k,  p) - ank_Mat52(n + 2 * k + 2, p) );

elseif strcmp(mode, 'Mat32')

    Ank = 1/2 * (ank_Mat32(n + 2 * k,  p) - ank_Mat32(n + 2 * k + 2, p) );

elseif contains(mode, 'Exp')

    if strcmp(mode, 'ExpNorm')

        Ank = 1/2 * (ank_ExpNorm(n + 2 * k,  p) - ank_ExpNorm(n + 2 * k + 2, p) );

    elseif contains(mode, 'ExpNonneg')

        Ank = 1/2 * (ank_ExpNonneg(n + 2 * k,  p) - ank_ExpNonneg(n + 2 * k + 2, p) );

    else

        Ank = 1/2 * (ank_Exp(n + 2 * k,  p) - ank_Exp(n + 2 * k + 2, p) );

    end

elseif strcmp(mode, 'SincNum')

    Ank = 1/2 * (ank_SincNum(n + 2 * k,  p) - ank_SincNum(n + 2 * k + 2, p) );

else

    Ank = 0;
    error('Unsupported mode');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ank] = ank_SqExp(n, p)

z = p^-2;
ank = pi * besseli(n, z, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ank] = ank_SqExpNorm(n, p)

norm_fac = 1 / ( 2 * pi * p^2 * (1 - exp(-2 / (p^2))));
ank = norm_fac * ank_SqExp(n, p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ank] = ank_SqExpNonneg(n, p)

ank = ank_SqExp(n, sqrt(2) * p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ank] = ank_Mat52(n, p)

z = 5 * (p^-2);
ank = (-1)^(n) *  pi * (hypergeomRegu([1], [1-n,  1+n], z)  -2 * z * hypergeomRegu(2, [2-n, 2+n], z) + 2 * z * sqrt(5)/p *  hypergeomRegu(2, [5/2-n, 5/2+n], z) + 2 * z * sqrt(pi) / 3 * hypergeomRegu([3/2, 2], [1/2, 2-n, 2+n], z) - z * sqrt(5*pi) / p *   hypergeomRegu([2, 5/2], [3/2, 5/2-n, 5/2+n], z)  );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ank] = ank_Mat32(n, p)

z = 3 * (p^-2);
ank = (-1)^(n) *  pi * (hypergeomRegu([1], [1-n,  1+n], z)  -2 * z * hypergeomRegu(2, [2-n, 2+n], z) + 2 * z * sqrt(3) / p *  hypergeomRegu(2, [5/2-n, 5/2+n], z)  );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ank] = ank_Exp(n, p)

ank = (-1)^(n) *  pi/(p) * (p * hypergeomRegu([1], [1-n,  1+n], p^-2)  - hypergeomRegu(1, [3/2-n, 3/2+n], p^-2) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ank] = ank_ExpNorm(n, p)

norm_fac = 1 / (2 * pi * (p^2 - p * exp(-2 / p) * (p + 2)));
ank = norm_fac * ank_Exp(n, p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ank] = ank_ExpNonneg(n, p)

ank = ank_Exp(n, 2 * p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ank] = ank_SincNum(n, p)

fun = @(x, n, p) cos(n * x) .* sinc(2 * sin(x / 2) * p / pi);
ank = integral(@(x) fun(x, n, p), 0, pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = hypergeomRegu(a, b, z)

y = 0;
for k = 0:100
%     y = y + exp( sum(logPorch(a, k)) + k * log(z) - (sum(log(1:k)) + sum(log(gamma(k+b))) ) );
%     y = y + real(exp( sum(logPorch(a, k)) + k * log(z) - (sum(log(1:k)) + sum(log(gamma(k+b))) ) ));
     y = y + real(exp( sum(logPorch(a, k)) + k * log(z) - (sum(log(1:k)) + sum(log( complex(gamma(k+b)) )) ) )); %Codegen
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = logPorch(x, n)

y = gammaln(x+n) - gammaln(x);

