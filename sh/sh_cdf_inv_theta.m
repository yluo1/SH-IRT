function [theta, stats] = sh_cdf_inv_theta(C_pdf, u, options)
%Compute inverse cumulative distribution function of density SH expansion over theta:
%theta = F_inv(u) where F(theta) is marginal CDF_{\Theta \Phi}(\theta, 2*pi) integrated over all of azimuth

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:                 [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%u:                 [N x 1] CDF targets

%options:           struct

%options.mode:      String, search method {'bisection', 'NewtonRaphsonBisection'}
%options.max_iter:  Maximum number of iterations
%options.tol:       zero-tolerance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%theta:             [N x M]

%stats:                   struct, performance stats
%stats.num_iter:          [N x M] Number of iterations
%stats.num_func_eval:     [N x M] Number of function evaluations
%stats.num_deriv_eval:    [N x M] Number of derivative evaluations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare forward and inverse CDFs

% rng(134);
% max_odr = 3;
% is_real = false; 
% C_pdf = sh_nrm(sh_msq(sh_rand(max_odr, 1, is_real), is_real), 'Sum'); %Is PDF

% N_fwd = 100;
% theta = linspace(0, pi, N_fwd)';
% N_inv = 11;
% CDF_C = sh_cdf_theta(C_pdf, theta);
% u = linspace(0, 1, N_inv)';
% theta_CDF_inv = sh_cdf_inv_theta(C_pdf, u);

% fontsize = 18;
% h_cdf = figure;
% h_cdf.Position = [100, 100, [600, 350] * 4/5];
% plot(theta, CDF_C, '-', theta_CDF_inv, u, '*', 'linewidth', 3, 'markersize', 10); grid on; axis tight; 
% h_lg = legend('$u = F_{\Theta}(\theta)$', '$\theta = F_{\Theta}^{-1}(u \sim \mathcal{U}(0,1))$', 'location', 'best', 'interpreter', 'latex');
% set(h_lg, 'fontsize', fontsize );
% xlabel('Co-latitude \theta', 'fontsize', fontsize);
% ylabel('CDF(\theta)', 'fontsize', fontsize);
% set(gca, 'fontsize', fontsize - 1);
% title('Cumulative Distribution Function', 'fontsize', fontsize + 1);

arguments
     C_pdf (:,:) double {coder.mustBeComplex} = complex(0);     
     u (:,1) double {mustBeNonnegative} = 0;

     options.mode (1,:) char {mustBeMember( options.mode, {'bisection', 'NewtonRaphsonBisection'} )} = 'NewtonRaphsonBisection';
     options.max_iter (1,1) double {mustBePositive} = 1000;
     options.tol (1,1) double {mustBePositive} = 1e-10;
end

assert(all(u >= 0) && all(u <= 1), 'u must be between [0, 1]');

[P, M] = size(C_pdf);
P = sqrt(P) - 1;
assert(P == floor(P), 'Invalid size C');

N = numel(u);

tol = options.tol;
max_iter = options.max_iter;

theta = zeros([N, M]);

if nargout > 1
    stats = struct('num_iter', zeros([N, M]), 'num_func_eval', zeros([N, M]), 'num_deriv_eval', zeros([N, M]));
end

for m = 1:M %Iterate over functions
    for n = 1:N %Iterate over CDF targets

        if strcmp(options.mode, 'bisection')

            up = pi;
            lo = 0;      
            if abs(sh_cdf_theta(C_pdf(:, m), lo) - u(n)) <= tol

                theta(n, m) = lo;

                if nargout > 1 %Track statistics
                     stats.num_iter(n, m) = 0;
                     stats.num_func_eval(n, m) = 1;
                end

            elseif  abs(sh_cdf_theta(C_pdf(:, m), up) - u(n)) <= tol

                theta(n, m) = up;

                if nargout > 1 %Track statistics
                     stats.num_iter(n, m) = 0;
                     stats.num_func_eval(n, m) = 2;
                end

            else

                theta(n, m) = (up + lo) / 2;
                for iter = 1:max_iter % Binary search
                    v = sh_cdf_theta(C_pdf(:, m), theta(n, m));
                    if abs(v - u(n)) <= tol
                        break;
                    elseif v - u(n) > 0
                        up = theta(n, m);
                    else
                        lo = theta(n, m);
                    end
                    theta(n, m) = (up + lo) / 2;
                end

                if nargout > 1 %Track statistics
                    stats.num_iter(n, m) = iter;
                    stats.num_func_eval(n, m) = 2 + iter;
                end
            end

        elseif strcmp(options.mode, 'NewtonRaphsonBisection')

            up = pi;
            lo = 0;      
            if abs(sh_cdf_theta(C_pdf(:, m), lo) - u(n)) <= tol

                theta(n, m) = lo;

                if nargout > 1 %Track statistics
                     stats.num_iter(n, m) = 0;
                     stats.num_func_eval(n, m) = 1;
                end

            elseif  abs(sh_cdf_theta(C_pdf(:, m), up) - u(n)) <= tol

                theta(n, m) = up;

                if nargout > 1 %Track statistics
                     stats.num_iter(n, m) = 0;
                     stats.num_func_eval(n, m) = 2;
                end

            else

                theta(n, m) = (up + lo) / 2;
                for iter = 1:max_iter % Newton-Raphson or binary search

                    [f_theta, d_f_theta] = sh_cdf_theta(C_pdf(:, m), theta(n, m));                                        
                    f_theta = f_theta - u(n);

                    if abs(f_theta) <= tol
                        break;
                    else

                        %Tighten bounds
                        if f_theta > 0
                            up = theta(n, m);
                        else
                            lo = theta(n, m);
                        end

                        %d_f_theta(1)

                        theta_NR = theta(n, m) - f_theta / d_f_theta; %Newton-Raphson theta
                        if theta_NR > lo && theta_NR < up %theta_NR within bounds

                            theta(n, m) = theta_NR;

                        else %theta_NR outside bounds, switch to bisection
                            theta(n, m) = (up + lo) / 2;
                        end

                    end
                end

                if nargout > 1 %Track statistics
                    stats.num_iter(n, m) = iter;
                    stats.num_func_eval(n, m) = 2 + iter;
                    stats.num_deriv_eval(n, m) = iter;
                end

                
            end

        else
            error('Unknown mode');
        end
    end
end