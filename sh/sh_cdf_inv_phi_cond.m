function [phi, stats] = sh_cdf_inv_phi_cond(C_pdf, u, theta, options)
%Compute inverse conditional cumulative distribution function 
%of spherical harmonic expansion probability density function
%of azimuth phi given co-latitude theta

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:                 [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%u:                 [N x 1] CDF targets
%theta:             [N x 1] Co-latitude between [0, pi]

%options:           struct

%options.mode:      String, search method {'bisection', 'NewtonRaphsonBisection'}
%options.max_iter:  Maximum number of iterations
%options.tol:       zero-tolerance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%phi:             [N x M]

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
% theta = pi / 3;
% phi = linspace(0, 2 * pi, N_fwd)';
% N_inv = 11;
% CDF_C = sh_cdf_phi_cond(C_pdf, theta, phi);
% u = linspace(0, 1, N_inv)';
% phi_CDF_inv = sh_cdf_inv_phi_cond(C_pdf, u, theta * ones(size(u)));

% fontsize = 18;
% h_cdf = figure;
% h_cdf.Position = [100, 100, [600, 350] * 4/5];
% plot(phi, CDF_C, '-', phi_CDF_inv, u, '*', 'linewidth', 3, 'markersize', 10); grid on; axis tight; 
% h_lg = legend('$u = F_{\Phi | \theta}(\phi)$', '$\phi = F_{\Phi | \theta}^{-1}(u \sim \mathcal{U}(0,1))$', 'location', 'best', 'interpreter', 'latex');
% set(h_lg, 'fontsize', fontsize );
% xlabel('Azimuth \phi', 'fontsize', fontsize);
% ylabel('CDF(\phi)', 'fontsize', fontsize);
% set(gca, 'fontsize', fontsize - 1);
% title('Cumulative Distribution Function', 'fontsize', fontsize + 1);

arguments
     C_pdf (:,:) double {coder.mustBeComplex} = complex(0);     
     u (:,1) double {mustBeNonnegative} = 0;
     theta (:,1) double {mustBeNonnegative} = pi/2;

     options.mode (1,:) char {mustBeMember( options.mode, {'bisection', 'NewtonRaphsonBisection'} )} = 'bisection';
     options.max_iter (1,1) double {mustBePositive} = 1000;
     options.tol (1,1) double {mustBePositive} = 1e-10;
end

assert(all(u >= 0) && all(u <= 1), 'u must be between [0, 1]');
assert(all(theta >= 0) && all(theta <= pi), 'theta must be between [0, pi]');
assert(numel(theta) == numel(u), 'theta, u size mismatch');

[P, M] = size(C_pdf);
P = sqrt(P) - 1;
assert(P == floor(P), 'Invalid size C');

N = numel(u);

tol = options.tol;
max_iter = options.max_iter;

phi = zeros([N, M]);

if nargout > 1
    stats = struct('num_iter', zeros([N, M]), 'num_func_eval', zeros([N, M]), 'num_deriv_eval', zeros([N, M]));
end

for m = 1:M %Iterate over functions
    for n = 1:N %Iterate over CDF targets

        if strcmp(options.mode, 'bisection')

            up = 2*pi;
            lo = 0;      
            if abs(sh_cdf_phi_cond(C_pdf(:, m), theta(n), lo) - u(n)) <= tol

                phi(n, m) = lo;

                if nargout > 1 %Track statistics
                     stats.num_iter(n, m) = 0;
                     stats.num_func_eval(n, m) = 1;
                end

            elseif  abs(sh_cdf_phi_cond(C_pdf(:, m), theta(n), up) - u(n)) <= tol

                phi(n, m) = up;

                if nargout > 1 %Track statistics
                     stats.num_iter(n, m) = 0;
                     stats.num_func_eval(n, m) = 2;
                end

            else

                phi(n, m) = (up + lo) / 2;
                for iter = 1:max_iter % Binary search
                    v = sh_cdf_phi_cond(C_pdf(:, m), theta(n), phi(n, m));
                    if abs(v - u(n)) <= tol
                        break;
                    elseif v - u(n) > 0
                        up = phi(n, m);
                    else
                        lo = phi(n, m);
                    end
                    phi(n, m) = (up + lo) / 2;
                end

                if nargout > 1 %Track statistics
                    stats.num_iter(n, m) = iter;
                    stats.num_func_eval(n, m) = 2 + iter;
                end
            end

        elseif strcmp(options.mode, 'NewtonRaphsonBisection')

            up = 2*pi;
            lo = 0;      
            if abs(sh_cdf_phi_cond(C_pdf(:, m), theta(n), lo) - u(n)) <= tol

                phi(n, m) = lo;

                if nargout > 1 %Track statistics
                     stats.num_iter(n, m) = 0;
                     stats.num_func_eval(n, m) = 1;
                end

            elseif  abs(sh_cdf_phi_cond(C_pdf(:, m), theta(n), up) - u(n)) <= tol

                phi(n, m) = up;

                if nargout > 1 %Track statistics
                     stats.num_iter(n, m) = 0;
                     stats.num_func_eval(n, m) = 2;
                end

            else

                phi(n, m) = (up + lo) / 2;
                for iter = 1:max_iter % Newton-Raphson or binary search

                    [f_phi, d_f_phi] = sh_cdf_phi_cond(C_pdf(:, m), theta(n), phi(n, m));                                        
                    f_phi = f_phi - u(n);

                    if abs(f_phi) <= tol
                        break;
                    else

                        %Tighten bounds
                        if f_phi > 0
                            up = phi(n, m);
                        else
                            lo = phi(n, m);
                        end

                        %d_f_phi(1)

                        phi_NR = phi(n, m) - f_phi / d_f_phi; %Newton-Raphson theta
                        if phi_NR > lo && phi_NR < up %theta_NR within bounds

                            phi(n, m) = phi_NR;

                        else %theta_NR outside bounds, switch to bisection
                            phi(n, m) = (up + lo) / 2;
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