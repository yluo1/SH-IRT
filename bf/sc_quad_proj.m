function [x, S, h_fig] = sc_quad_proj(a, b, mode, options)
% Find zeros of the secular quadratic projection equation:
% S(x) = sum_n ( a(n) * b(n) ) / ( x - b(n) )^2 = 0 for squared modulus a = u .* conj(u), real b
% Critical zero is real zero x that satisfies  b(b < 0) < x < b(b > 0)

% mode = 'polysym'                   Convert secular equation to symbolic univariate polynomial and find all zeros
% mode = 'bisection'                 Bisection method to find critical zero within bounds
% mode = 'NewtonRaphson'             Newton-Raphson and bisection method when out-of-bounds to find critical zero

% Paper reference:
% Luo, Y. (2024). Constant Directivity Loudspeaker Beamforming. In 2024 32nd European Signal Processing Conference (EUSIPCO) (pp. 246-250). IEEE.

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%a:         [1 x N] Non-negative
%b:         [1 x N] Real
%mode:      String, compute method {'polysym', 'bisection', 'NewtonRaphsonBisection'}

%options:               struct
%options.enable_disp:   Logical, if true, plot secular equation and zeros

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x:         [1 x 2 * (N - 1)] zeros if mode = 'polysym', otherwise [1 x 1] critical zero
%S:         [1 x 2 * (N - 1)] or [1 x 1] S(x)
%h_fig:     Figure handle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compute zeros and plot secular equation across modes

%rng(2342);
%u = [12, 1, 10];
%a = u .* conj(u);
%b = [-2, -1,  1];
%[x_polysym, S_polysym]     = sc_quad_proj(a, b, 'polysym', 'enable_disp', true); x_polysym
%[x_bisection, S_bisection] = sc_quad_proj(a, b, 'bisection', 'enable_disp', true); x_bisection
%[x_nr, S_nr] = sc_quad_proj(a, b, 'NewtonRaphson', 'enable_disp', true); x_nr

arguments
    a (1,:) double {mustBeNonnegative} = 1;
    b (1,:) double {mustBeReal} = 1;
    mode (1,:) char {mustBeMember(mode, {'polysym', 'bisection', 'NewtonRaphson'})} = 'NewtonRaphson';

    options.enable_disp (1,1) logical = false;
    options.display_deriv (1,1) logical = false;
    options.display_pos_neg (1,1) logical = false;
    options.disp_critical (1,1) logical = true;
    options.disp_position (1,4) double = [100, 100, 600, 400];
    options.fontsize (1,1) {mustBePositive} = 20;
end

N = numel(a);
assert(N == numel(b), 'a, b size mismatch');

S_func = @(x) sum((a .* b ) ./ (x(:) - b).^2, 2);              % [1 x numel(x)] secular equation
dSdx_func  = @(x) -2 * sum(( a .* b ) ./ (x(:) - b).^3, 2);    % [1 x numel(x)] first derivative
dSdx2_func = @(x) 6  * sum(( a .* b ) ./ (x(:) - b).^4, 2);    % [1 x numel(x)] second derivative

if strcmp(mode, 'polysym')  % Symbolic univariate polynomial root-finding

    syms lambda
    p = sym(0);
    for n = 1:N
        p_n = sym( a(n) * b(n) );
        for i = 1:N
            if i ~= n            
                p_n = p_n * sym( (lambda -  b(i))^2 );
            end
        end
        p = p + p_n;
    end    
    x = double(root(p, lambda));

elseif strcmp(mode, 'bisection') % Bisection method

    [x, exitflag, stats] = bisect_zero(S_func, max(b(b<0)), min(b(b>0)), 'bisection');
    % exitflag
    % stats

elseif strcmp(mode, 'NewtonRaphson') % Newton Rhapson with Bisection when out-of-bounds

    [x, exitflag, stats] = bisect_zero(S_func, max(b(b<0)), min(b(b>0)), 'NewtonRaphson', 'd_func', dSdx_func);
    % exitflag
    % stats

else
    error('Unknown mode');
end
S = S_func(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_fig = [];
if options.enable_disp && coder.target("MATLAB")
    
    % Ascended sorted b to separate into pos/neg halves
    [b_asc, idx] = sort(b, 'ascend');
    a_asc = a(idx);

    mask_pos = b_asc > 0;
    mask_neg = b_asc < 0;

    idx_pos = find(mask_pos, 1, 'first'); % Index first positive value coeffs b
    idx_neg = idx_pos - 1;
    
    x_real_zeros = x(imag(x) == 0);
    idx_crit = find(x_real_zeros > b_asc(idx_neg) & x_real_zeros < b_asc(idx_pos), 1, 'first'); % Critical zero

    x_disp = linspace(min([b(:); x_real_zeros(:)]) - 1, max([b(:); x_real_zeros(:)]) + 1, 1024)';
    y_disp = S_func(x_disp);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h_fig = figure;
    h_fig.Position = options.disp_position;
  
    yline(0, 'k--', 'HandleVisibility', 'off');  hold on;
    
    lg_str = {};
    plot(x_disp, y_disp, 'b-', 'linewidth', 2);                lg_str{end+1} = '$S(\lambda)$';

    if options.display_deriv    % Display derivatives

        y_disp_deriv  = dSdx_func(x_disp);
        y_disp_deriv2 = dSdx2_func(x_disp);
        plot(x_disp, y_disp_deriv, 'r-.', 'linewidth', 1.5);   lg_str{end+1} = '$\frac{d S(\lambda)}{d \lambda}$';
        plot(x_disp, y_disp_deriv2, 'g-', 'linewidth', 1.5);   lg_str{end+1} = '$\frac{d^2 S(\lambda)}{d \lambda}$';
    end
    
    if options.display_pos_neg  % Display positive and negative function components

        S_pos_func = @(x) sum((a_asc(mask_pos) .* b_asc(mask_pos) ) ./ (x(:) - b_asc(mask_pos) ).^2, 2);
        S_neg_func = @(x) sum((a_asc(mask_neg) .* b_asc(mask_neg) ) ./ (x(:) - b_asc(mask_neg) ).^2, 2);
        y_pos_plot = S_pos_func(x_disp);
        y_neg_plot = S_neg_func(x_disp);
        
        plot(x_disp, y_neg_plot, '--', 'color', [0.9294, 0.6941, 0.1294], 'linewidth', 2);     lg_str{end+1} = '$S_{-}(\lambda)$';
        plot(x_disp, y_pos_plot, '--', 'color', [0.4941, 0.1843, 0.5569], 'linewidth', 2);     lg_str{end+1} = '$S_{+}(\lambda)$';
    end

    % Display poles and zeros
    for n = 1:numel(b)  % Poles
        if n == 1
            plot(b(n), 0, '*', 'color', [0, 0.4471, 0.7412], 'linewidth', 2, 'markersize', 10);
            lg_str{end+1} =  'Poles';
        else
            plot(b(n), 0, '*', 'color', [0, 0.4471, 0.7412], 'linewidth', 2, 'markersize', 10, 'HandleVisibility', 'off');
        end
    end    
    for n = 1:numel(x_real_zeros)  % Zeros      
        if n == 1
            plot(x_real_zeros(n), 0, 'ms', 'linewidth', 2, 'markersize', 12);
            lg_str{end+1} =  'Zeros';
        else
            plot(x_real_zeros(n), 0, 'ms', 'linewidth', 2, 'markersize', 12, 'HandleVisibility', 'off');
        end
    end

    if options.disp_critical % Display critical zero

        plot(x_real_zeros(idx_crit), 0, 'ro', 'linewidth', 1.65, 'markersize', 16);                                 lg_str{end+1} = '$\lambda_*$';

        xline(b_asc(idx_neg), '-.', 'color', [0.9294, 0.6941, 0.1294], 'linewidth', 2);                             lg_str{end+1} = '$b_-$';
        xline(2*b_asc(idx_neg) - x_real_zeros(idx_crit), ':', 'color', [0.9294, 0.6941, 0.1294], 'linewidth', 2);   lg_str{end+1} = '$2 b_- - \lambda_*$';

        xline(b_asc(idx_pos), '-.', 'color', [0.4941, 0.1843, 0.5569], 'linewidth', 2);                             lg_str{end+1} = '$b_+$';
        xline(2*b_asc(idx_pos) - x_real_zeros(idx_crit), ':', 'color', [0.4941, 0.1843, 0.5569], 'linewidth', 2);   lg_str{end+1} = '$2 b_+ - \lambda_*$';
        
    end
    
    grid on;
    xlim( [min(x_disp), max(x_disp)] );
    prc = prctile( y_disp, [15, 85] ); % Compute percentile
    ylim( [min(prc(1), -10), max(prc(2), 10)] );

    xlabel('\lambda', 'fontsize', options.fontsize);
    ylabel('S(\lambda)', 'fontsize', options.fontsize);
    title('Secular Quadratic Projection Equation', 'fontsize', options.fontsize + 1);
    set(gca, 'fontsize', options.fontsize - 1);
    lg = legend(lg_str, 'location', 'best', 'interpreter', 'latex', 'NumColumns', 2);
    set(lg, 'fontsize', options.fontsize);

end