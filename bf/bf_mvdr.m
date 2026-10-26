function [W, fval, WNG] = bf_mvdr(v, R, kappa, mode)
% Minimum variance distortionless response beamformer design with condition number constraints
% min_{W, sigma} W'*R(kappa)*W   s.t.    W'*v = 1,   

% mode = 'diag'   where R is diagonal loaded R(kappa) = R + sigma * I  s.t. cond(R + sigma * I) <= kappa 
% =>    (lambda_max + sigma) / (lambda_min + sigma) <= kappa
% =>    sigma >= (lambda_max - lambda_min) / (kappa - 1)
% =>    sigma = (lambda_max - lambda_min) / (kappa - 1) at minimizer

% mode = 'thres' where  R(kappa) = V*D*V' for eigenvector matrix V s.t. cond(D) <= kappa 
% Eigenvalues below threshold are limited:
% D = diag(max(lambda, lambda_max/kappa))

% mode = 'exp'  where  R(kappa) = V*D*V' for eigenvector matrix V s.t. cond(D) <= kappa 
% Eigenvalues are exponentiated s.t. cond(D) <= kappa 
% D = diag(lambda.^p)
% Solves for max_p  s.t.  (lambda_max / lambda_min)^p <= kappa,      0 <= p <= 1
% =>    p = min(1, log(kappa) / (log(lambda_max) - log(lambda_min)) )

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%v:         [N x 1]     Steering vector for N channels
%R:         [N X N]     Rejection covariance matrix
%kappa:                 Scalar, condition number constraint >= 1, inf = unconstrained
%                       Choose kappa <= SNR
%mode:                  String, constraint method {'diag', 'thres', 'trunc'}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%W:         [N x 1]     Conjugate symmetric beamformer weights
%fval                   Scalar, objective real(W'*R*W)
%WNG:                   Scalar, white-noise gain W'*W

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Plot condition number versus WNG tradeoff for poorly conditioned R

% rng(4113);
% N = 6;
% v = randn(N, 1) + randn(N, 1)*1i;;
% R = randn(N) + randn(N)*1i; R = R*R';
% R = R - 0.999 * eigs(R, 1, 'smallestreal') * eye(N);
% N_kappa = 32;
% kappa_list = logspace(log10(1.1), log10(1000000), N_kappa);
% W_diag = zeros(N, N_kappa);    WNG_diag = zeros(1, N_kappa);
% W_thres = zeros(N, N_kappa);   WNG_thres = zeros(1, N_kappa);
% W_exp = zeros(N, N_kappa);     WNG_exp = zeros(1, N_kappa);
% for n = 1:N_kappa
%   [W_diag(:, n), ~, WNG_diag(n)]     = bf_mvdr(v, R, kappa_list(n), 'diag');
%   [W_thres(:, n), ~, WNG_thres(n)]   = bf_mvdr(v, R, kappa_list(n), 'thres');
%   [W_exp(:, n), ~, WNG_exp(n)]       = bf_mvdr(v, R, kappa_list(n), 'exp');
% end
% 
% %Plotting
% fontsize = 14;
% figure;
% semilogx(kappa_list, mag2db(WNG_diag), kappa_list, mag2db(WNG_thres), kappa_list, mag2db(WNG_exp), 'linewidth', 1.5); 
% grid on; axis tight;
% xlabel('Max Condition Number \kappa Constraint', 'fontsize', fontsize);
% ylabel('White-noise Gain (dB)', 'fontsize', fontsize);
% title('MVDR White-Noise Gain across Conditioning Methods', 'fontsize', fontsize + 1); set(gca, 'fontsize', fontsize - 1);
% h_lg = legend('diag', 'thres', 'exp', 'location', 'best'); set(h_lg, 'fontsize', fontsize - 1);

arguments
    v (:,1) double {coder.mustBeComplex} = complex(1);
    R (:,:) double {coder.mustBeComplex} = complex(1);
    kappa (1,1) double {mustBePositive} = inf;
    mode (1,:) char {mustBeMember(mode, {'diag', 'thres', 'exp'} )} = 'diag';
end

N = numel(v);

assert(N == size(R, 1) && N == size(R, 2), 'v and R size mismatch');
assert(kappa >= 1, 'Require condition number kappa >= 1');

if strcmp(mode, 'diag') % Shape preserving across entire eigenvalues by upwards shifting all eigenvalues
    % minimizes total difference of normalized cumulative eigenvalues and constrained eigenvalues

    lambda_max = eigs(R, 1, 'largestreal');
    lambda_min = eigs(R, 1, 'smallestreal');

    sigma = (lambda_max - lambda_min) / (kappa - 1);
    Rinv_v = (R + sigma * eye(N)) \ v; 

elseif strcmp(mode, 'thres') % Shape preserving upto threshold, then brick-wall limit small eigenvalues
    % minimizes total difference between original and constrained eigenvalues

    [V, D] = eig(R);
    d = diag(D);
    [lambda, idx] = sort(d, 'descend');
    V = V(:, idx);
    
    lambda_max = lambda(1);
    Dinv = diag( 1./ max(lambda, lambda_max/kappa) );
    
    Rinv_v = (V * Dinv * V') * v;


elseif strcmp(mode, 'exp') % Shape smoothing eigenvalues via soft-compression

    [V, D] = eig(R);
    d = diag(D);
    [lambda, idx] = sort(d, 'descend');
    V = V(:, idx);
    
    lambda_max = lambda(1);
    lambda_min = lambda(end);
    p = min(1, log(kappa) / (log(lambda_max) - log(lambda_min)));

    Dinv = diag( lambda.^(-p) );
    
    Rinv_v = (V * Dinv * V') * v;

else
    error('Unknown mode');
end

% Output
W = Rinv_v / (v' * Rinv_v); % Normalize for distortionless constraint
fval = real(W'*R*W);
WNG = W'*W;
