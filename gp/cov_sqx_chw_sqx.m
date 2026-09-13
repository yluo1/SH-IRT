function [K, dK_dsigma, dK_dell_c, dK_dell_f] = cov_sqx_chw_sqx(X, Y, sigma, ell_c, ell_f)
%Compute square exponential of chordal distance x squared exponential of log-frequency distance

%K = sigma^2 * exp(-d(theta_x, phi_x, theta_y, phi_y)^2 / ( 2 * ell_c^2 ))
%            * exp(-(log(f_x), log(f_y))^2 / (2 * ell_f^2) )^2 )

%theta:      Co-latitude [0, pi]
%phi:        Azimuth [0, 2 * pi)
%omega:      Angular frequency
%f:          Ordinary frequency

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X:        [NX x 3]  Row-matrix of [f_x, theta_x, phi_x]
%Y:        [NY x 3]  Row-matrix of [f_y, theta_y, phi_y]

%sigma:    Scalar, covariance scaling hyperparameter
%ell_c:    Scalar, chordal distance scale hyperparameter
%ell_f:    Scalar, log-frequency scale hyperparameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%K:             [NX x NY] Covariance matrix

%dK_dsigma:     [NX x NY] Partial derivative matrix of entrants in K w.r.t. sigma
%dK_dell_c:     [NX x NY] Partial derivative matrix of entrants in K w.r.t. ell_c
%dK_dell_f:     [NX x NY] Partial derivative matrix of entrants in K w.r.t. ell_f

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Evaluate over random grid

% rng(521);
% 
% N = 20;
% X = [2 * pi * 343 ./ rand(N,1), rand(N,1) * pi, rand(N,1) * 2*pi];
% sigma = 1;
% ell = 0.3;
% gamma = 0;
% 
% K = cov_sqx_chw_sqx(X, X, sigma, ell, gamma);
% eig_K = eig((K+K')/2)
% min(eig_K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Evaluate over uniform grid

% N = 10;
% [theta, phi] = sh_fib(N);
% X = [2 * pi * ones(N,1) * 200, theta, phi];
% sigma = 1;
% ell = 0.3;
% gamma = 0;
% 
% K = cov_sqx_chw_sqx(X, X, sigma, ell, gamma);
% eig_K = eig((K+K')/2)
% min(eig_K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Evaluate over identical spherical coordinates, varying frequency

% N = 10;
% X = [2 * pi * logspace(log10(50), log10(24000), N)', ones(N, 1) * 0.35, ones(N, 1) * 1.24];
% sigma = 1;
% ell = 0.3;
% gamma = 0;
% 
% K = cov_sqx_chw_sqx(X, X, sigma, ell, gamma);
% eig_K = eig((K+K')/2)
% min(eig_K)


arguments
    X (:,3) double = [100 0 0];
    Y (:,3) double = [100 0 0];
    sigma (1,1) double {mustBePositive} = 1;
    ell_c (1,1) double {mustBePositive} = 1;
    ell_f (1,1) double {mustBePositive} = 1;
end

NX = size(X, 1);
NY = size(Y, 1);

%Compute chordal distance
VX = zeros([NX, 3]);
[VX(:,1), VX(:,2), VX(:,3)] = sph2cart( X(:,3), pi/2 - X(:,2), ones(NX, 1) );

VY = zeros([NY, 3]);
[VY(:,1), VY(:,2), VY(:,3)] = sph2cart( Y(:,3), pi/2 - Y(:,2), ones(NY, 1) );

%Compute chordal distance
%D = 2 * sin(abs( acos(max(min(VX * VY', 1), -1) )  ) / 2);
%D = pdist2(VX, VY);
D_sq = pdist2(VX, VY, 'squaredeuclidean');

%Compute log-frequency distance
f_x = X(:,1) / (2 * pi);
f_y = Y(:,1) / (2 * pi);
F = repmat(log(f_x), [1, NY]) - repmat(log(f_y'), [NX, 1]);
F_sq = F.^2;

%Squared exponential term
xterm_c = exp( - D_sq ./ (2 * ell_c^2) );
xterm_f = exp( - F_sq ./ (2 * ell_f^2) );

K =  sigma^2 .* xterm_c .* xterm_f;

%Compute partial derivatives
if nargout > 1
    dK_dsigma = 2 * K / sigma;
end

if nargout > 2
    dK_dell_c = K .* D_sq ./ ell_c^3;
end

if nargout > 3
    dK_dell_f = K .* F_sq ./ ell_f^3;
end

