function [K, dK_dsigma, dK_dell] = cov_sqx_chw(X, Y, sigma, ell)
%Compute square exponential of chordal distance

%K =  sigma^2  *  exp(-d(theta_x, phi_x, theta_y, phi_y)^2 / ( 2 * ell^2 ))

%theta:      Co-latitude [0, pi]
%phi:        Azimuth [0, 2 * pi)
%sigma:      Standard deviation for covariance scaling

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X:        [NX x 2]  Row-matrix of [theta_x, phi_x]
%Y:        [NY x 2]  Row-matrix of [theta_y, phi_y]

%sigma:    Scalar, covariance scaling hyperparameter
%ell:      Scalar, wavelength scale hyperparameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%K:             [NX x NY] Covariance matrix

%dK_dsigma:     [NX x NY] Partial derivative matrix of entrants in K w.r.t. sigma
%dK_dell:       [NX x NY] Partial derivative matrix of entrants in K w.r.t. ell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Evaluate over random grid

% rng(521);
% 
% N = 20;
% X = [rand(N,1) * pi, rand(N,1) * 2*pi];
% sigma = 1;
% ell = 0.3;
% 
% K = cov_sqx_chw(X, X, sigma, ell);
% eig_K = eig((K+K')/2)
% min(eig_K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Evaluate over uniform grid

% N = 10;
% [theta, phi] = sh_fib(N);
% X = [theta, phi];
% sigma = 1;
% ell = 0.3;
% 
% K = cov_sqx_chw(X, X, sigma, ell);
% eig_K = eig((K+K')/2)
% min(eig_K)


arguments
    X (:,2) double = [0 0];
    Y (:,2) double = [0 0];
    sigma (1,1) double {mustBePositive} = 1;
    ell (1,1) double {mustBePositive} = 1;
end

NX = size(X, 1);
NY = size(Y, 1);

%Compute chordal distance
VX = zeros([NX, 3]);
[VX(:,1), VX(:,2), VX(:,3)] = sph2cart( X(:,2), pi/2 - X(:,1), ones(NX, 1) );

VY = zeros([NY, 3]);
[VY(:,1), VY(:,2), VY(:,3)] = sph2cart( Y(:,2), pi/2 - Y(:,1), ones(NY, 1) );

%Compute distance
%D = 2 * sin(abs( acos(max(min(VX * VY', 1), -1) )  ) / 2);
%D = pdist2(VX, VY);
%D_sq = D.^2;
D_sq = pdist2(VX, VY, 'squaredeuclidean');

%Squared exponential term
xterm = exp( - D_sq ./ (2 * ell^2) );

K =  sigma^2 .* xterm;

%Compute partial derivatives
if nargout > 1
    dK_dsigma = 2 * K / sigma;
end

if nargout > 2
    dK_dell = K .* D_sq ./ ell^3;
end
