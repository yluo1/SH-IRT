function [K, dK_dsigma, dK_dell, dK_dgamma] = cov_sqx_chw_ns(X, Y, sigma, ell, gamma)
%Compute squared exponential of chordal distance with non-stationary wavelength covariance matrix

%K =  sigma^2  * lambda_x^(1/2) * lambda_y^(1/2) * ((lambda_x^2 + lambda_y^2 ) / 2)^(-1/2) 
%  *  exp(-d(theta_x, phi_x, theta_y, phi_y)^2 / ( (lambda_x^2 + lambda_y^2) / 2) )

%theta:      Co-latitude [0, pi]
%phi:        Azimuth [0, 2 * pi)
%sigma:      Standard deviation for covariance scaling
%lambda:     Scaled wavelength (lambda = ell * unscaled_lambda^gamma) (0, inf)

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X:        [NX x 3]  Row-matrix of [unscaled_lambda_x, theta_x, phi_x]
%Y:        [NY x 3]  Row-matrix of [unscaled_lambda_y, theta_y, phi_y]

%sigma:    Scalar, covariance scaling hyperparameter
%ell:      Scalar, velocity hyperparameter
%gamma:    Scalar, frequency power hyperparameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%K:             [NX x NY] Covariance matrix

%dK_dsigma:     [NX x NY] Partial derivative matrix of entrants in K w.r.t. sigma
%dK_dell:       [NX x NY] Partial derivative matrix of entrants in K w.r.t. ell
%dK_dgamma:     [NX x NY] Partial derivative matrix of entrants in K w.r.t. gamma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Evaluate over random grid

% rng(521);
% 
% N = 20;
% X = [rand(N,1) * 343, rand(N,1) * pi, rand(N,1) * 2*pi];
% sigma = 1;
% ell = 0.3;
% gamma = 0;
% 
% K = cov_sqx_chw_ns(X, X, sigma, ell, gamma);
% eig_K = eig((K+K')/2)
% min(eig_K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Evaluate over uniform grid

% N = 10;
% [theta, phi] = sh_fib(N);
% X = [ones(N,1) * 10, theta, phi];
% sigma = 1;
% ell = 0.3;
% gamma = 0;
% 
% K = cov_sqx_chw_ns(X, X, sigma, ell, gamma);
% eig_K = eig((K+K')/2)
% min(eig_K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Evaluate over identical spherical coordinates, varying wavelengths

% N = 10;
% X = [(1:N)' * 10, ones(N, 1) * 0.35, ones(N, 1) * 1.24];
% sigma = 1;
% ell = 0.3;
% gamma = 0;
% 
% K = cov_sqx_chw_ns(X, X, sigma, ell, gamma);
% eig_K = eig((K+K')/2)
% min(eig_K)


arguments
    X (:,3) double = [1 0 0];
    Y (:,3) double = [1 0 0];
    sigma (1,1) double {mustBePositive} = 1;
    ell (1,1) double {mustBePositive} = 1;
    gamma (1,1) double {mustBeNonnegative} = 1;
end

NX = size(X, 1);
NY = size(Y, 1);

%Compute chordal distance
VX = zeros([NX, 3]);
[VX(:,1), VX(:,2), VX(:,3)] = sph2cart( X(:,3), pi/2 - X(:,2), ones(NX, 1) );

VY = zeros([NY, 3]);
[VY(:,1), VY(:,2), VY(:,3)] = sph2cart( Y(:,3), pi/2 - Y(:,2), ones(NY, 1) );

%Compute distance
D = 2 * sin(abs( acos(max(min(VX * VY', 1), -1) )  ) / 2);

%Unscaled wavelengths 
unscaled_lambda_x = X(:, 1);
unscaled_lambda_y = Y(:, 1);
unscaled_lambda_x_mat = repmat(unscaled_lambda_x,  [1, NY]);
unscaled_lambda_y_mat = repmat(unscaled_lambda_y', [NX, 1]);

%Scaled wavelengths
X(:, 1) = ell * (X(:, 1).^gamma);
Y(:, 1) = ell * (Y(:, 1).^gamma);

%Compute non-stationary wavelength length-scales
L = (repmat(X(:, 1).^2, [1, NY]) + repmat((Y(:, 1).^2)', [NX, 1])) / 2;

%Normalization term
c = repmat( sqrt(X(:, 1)),  [1, NY]) ...
 .* repmat( sqrt(Y(:, 1))', [NX, 1]) ...
 ./ sqrt(L);                 

%Squared exponential term
xterm = exp( - (D.^2) ./ L );

K =  sigma^2  * c .* xterm;

%Compute partial derivatives
if nargout > 1
    dK_dsigma = 2 * K / sigma;
end

if nargout > 2

    %K =  sigma^2  * lambda_x^(1/2) * lambda_y^(1/2) * ((lambda_x^2 + lambda_y^2 ) / 2)^(-1/2) 
    %  *  exp(-d(theta_x, phi_x, theta_y, phi_y)^2 / ( (lambda_x^2 + lambda_y^2) / 2) )
    
    % lambda = ell * unscaled_lambda^gamma
    % dlambda_x_dell = unscaled_lambda_x^gamma
    % dlambda_y_dell = unscaled_lambda_y^gamma
    
    %lnorm = sqrt(lambda_x * lambda_y / ((lambda_x^2 + lambda_y^2 ) / 2)) 
    %      = sqrt(2 * (unscaled_lambda_x * unscaled_lambda_y)^gamma / (unscaled_lambda_x^(2 * gamma) + unscaled_lambda_y^(2 * gamma)) )       
    
    f = (unscaled_lambda_x_mat .* unscaled_lambda_y_mat).^gamma;
    g = unscaled_lambda_x_mat.^(2 * gamma) + unscaled_lambda_y_mat.^(2 * gamma);

    %lnorm = sqrt(2 * (unscaled_lambda_x_mat .* unscaled_lambda_y_mat).^gamma ./ (unscaled_lambda_x_mat.^(2 * gamma) + unscaled_lambda_y_mat.^(2 * gamma)) );
    lnorm = sqrt(2 * f ./ g);

    %dxterm_dell = exp(-d(theta_x, phi_x, theta_y, phi_y)^2 / ( ell^2 * (unscaled_lambda_x^(2 * gamma) + unscaled_lambda_y^(2 * gamma)) / 2) )
    %            *  4 * d(theta_x, phi_x, theta_y, phi_y)^2 / ( ell^3 * (unscaled_lambda_x^(2 * gamma) + unscaled_lambda_y^(2 * gamma)))           
    
    dxterm_dell = 4 * xterm .* (D.^2) ./ (ell^3 * g);
    
    %K = sigma^2 .* lnorm .*  xterm
    dK_dell = sigma^2 * (lnorm .* dxterm_dell);

end

if nargout > 3

    %K =  sigma^2  * lambda_x^(1/2) * lambda_y^(1/2) * ((lambda_x^2 + lambda_y^2 ) / 2)^(-1/2) 
    %  *  exp(-d(theta_x, phi_x, theta_y, phi_y)^2 / ( (lambda_x^2 + lambda_y^2) / 2) )
    
    % lambda = ell * unscaled_lambda^gamma
    % dlambda_x_dgamma = ell * log(gamma) * unscaled_lambda_x^gamma
    % dlambda_y_dgamma = ell * log(gamma) * unscaled_lambda_y^gamma
    
    %lnorm = sqrt(lambda_x * lambda_y / ((lambda_x^2 + lambda_y^2 ) / 2)) 
    %      = sqrt(2 * (unscaled_lambda_x * unscaled_lambda_y)^gamma / (unscaled_lambda_x^(2 * gamma) + unscaled_lambda_y^(2 * gamma)) )       

      
    df_dgamma = log(unscaled_lambda_x_mat .* unscaled_lambda_y_mat) .* f;
    dg_dgamma = 2 * (log(unscaled_lambda_x_mat) .* unscaled_lambda_x_mat.^(2 * gamma) + log(unscaled_lambda_y_mat) .* unscaled_lambda_y_mat.^(2 * gamma));
    

    dlnorm_dgamma = (1./lnorm) .* (df_dgamma .* g - f .* dg_dgamma) ./ (g.^2);
           
    %dxterm_dgamma = exp(-d(theta_x, phi_x, theta_y, phi_y)^2 / ( ell^2 * (unscaled_lambda_x^(2 * gamma) + unscaled_lambda_y^(2 * gamma)) / 2) )
    %              *  2 / ell^2 * d(theta_x, phi_x, theta_y, phi_y)^2 / ( 2 * (log(unscaled_lambda_x) .* unscaled_lambda_x^(2 * gamma) + log(unscaled_lambda_y) .* unscaled_lambda_y^(2 * gamma)))           
    
    dxterm_dgamma = (2 / ell^2) .* xterm .* (D.^2) .* g.^(-2)  .* dg_dgamma;
    
    %K = sigma^2 .* lnorm .*  xterm
    dK_dgamma = sigma^2 * (lnorm .* dxterm_dgamma + dlnorm_dgamma .* xterm);

end

