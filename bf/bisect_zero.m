function [x, exitflag, stats] = bisect_zero(func, lb, ub, mode, options)
% Bisection method variants for root-finding func(x) in an interval [lb, ub] containing the zero
% Requirements: func(lb) <= 0 <= func(ub) or func(lb) >= 0 >= func(ub)

% mode = 'bisection':           Bisection method
% mode = 'NewtonRaphson':       Newton-Raphson when x is inside bound (requires options.d_func), 
%                               bisection method when x goes outside bound

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%func:                  Function f(x)
%lb:                    Lower bound of x
%ub:                    Upper bound of x
%mode:                  String, method {'bisection', 'NewtonRaphson'}

%options:               struct
%options.d_func:        Derivative function
%options.tol:           Scalar, zero tolerance for convergence
%options.max_iter:      Maximum number of iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%x:                     Scalar, zero
%exitflag:              1 = solution found within options.max_iter, 0 otherwise
%stats:                 struct containing number of function evaluations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample usage: Find root of cubic polynomial

%func = @(x) x.^3 - 2*x + 2;
%d_func = @(x) 3*x.^2 - 2;
%lb = -10;
%ub = 10;
%[x_bis, exitflag_bis, stats_bis] = bisect_zero(func, lb, ub, 'bisection')
%[x_NR, exitflag_NR, stats_NR] = bisect_zero(func, lb, ub, 'NewtonRaphson', 'd_func', d_func)
%func(x_bis)
%func(x_NR)

% %Negated function
%func = @(x) -(x.^3 - 2*x + 2);
%d_func = @(x) -(3*x.^2 - 2);
%lb = -10;
%ub = 10;
%[x_bis, exitflag_bis, stats_bis] = bisect_zero(func, lb, ub, 'bisection')
%[x_NR, exitflag_NR, stats_NR] = bisect_zero(func, lb, ub, 'NewtonRaphson', 'd_func', d_func)
%func(x_bis)
%func(x_NR)

arguments
    func function_handle {mustBeNonempty} = @(x) x
    lb (1,1) double {mustBeFinite} = 1;
    ub (1,1) double {mustBeFinite} = -1;
    mode (1,:) char {mustBeMember(mode, {'bisection', 'NewtonRaphson'} )} = 'NewtonRaphson';

    options.d_func  function_handle = function_handle.empty
    options.tol (1,1) double {mustBePositive} = 1e-10;
    options.max_iter (1,1) double {mustBePositive} = 1000;
end

assert(lb <= ub, 'lb must be <= ub');

tol = options.tol;
max_iter = options.max_iter;
exitflag = 0;

if nargout > 2
    stats = struct('num_iter', 0, 'num_func_eval', 0, 'num_deriv_eval', 0);
end

f_lb = func(lb);
f_ub = func(ub);
assert(sign(f_lb) * sign(f_ub) <= 0, 'func(lb) and func(ub) must bound 0');

if f_lb >= 0 && f_ub <= 0 % Negate the function such that func(lb) <= 0 <= func(ub)
    func = @(x) -func(x);
    if ~isempty(options.d_func)
        options.d_func = @(x) -options.d_func(x);
    end
    f_lb = -f_lb;
    f_ub = -f_ub;
end

% Following methods assume func(lb) <= 0 <= func(ub)
if strcmp(mode, 'bisection')
   
    if abs(f_lb) <= tol
    
        x = lb;
        if nargout > 2 % Track statistics
             stats.num_iter = 0;
             stats.num_func_eval = 2;
        end
        exitflag = 1;

    elseif abs(f_ub) <= tol
    
        x = ub;
        if nargout > 2 % Track statistics
             stats.num_iter = 0;
             stats.num_func_eval = 2;
        end
        exitflag = 1;

    else
    
        x = (ub + lb) / 2;
        for iter = 1:max_iter % Binary search
            v = func(x);
            if abs(v) <= tol % Converged
                exitflag = 1;
                break;
            elseif v > 0
                ub = x;
            else
                lb = x;
            end
            x = (ub + lb) / 2;
        end
    
        if nargout > 2 % Track statistics
            stats.num_iter = iter;
            stats.num_func_eval = 2 + iter;
        end
    end

elseif strcmp(mode, 'NewtonRaphson')
     
    assert(~isempty(options.d_func), 'NewtonRaphson requires options.d_func');

    if abs(f_lb) <= tol
    
        x = lb;
        if nargout > 2 % Track statistics
             stats.num_iter = 0;
             stats.num_func_eval = 2;
        end
        exitflag = 1;

    elseif abs(f_ub) <= tol
    
        x = ub;
        if nargout > 2 % Track statistics
             stats.num_iter = 0;
             stats.num_func_eval = 2;
        end
        exitflag = 1;
    
    else
    
        x = (ub + lb) / 2;
        for iter = 1:max_iter % Newton-Raphson or binary search
    
            f_x = func(x);
            df_dx = options.d_func(x);
    
            if abs(f_x) <= tol % Converged
                exitflag = 1;
                break;
            else
    
                % Tighten bounds
                if f_x > 0
                    ub = x;
                else
                    lb = x;
                end
    
                x_NR = x - f_x / df_dx;   % Newton-Raphson step
                if x_NR > lb && x_NR < ub % x_NR within bounds
                    x = x_NR;
                else % x_NR outside bounds, switch to bisection
                    x = (ub + lb) / 2;
                end
    
            end
        end
    
        if nargout > 2 %Track statistics
            stats.num_iter = iter;
            stats.num_func_eval = 2 + iter;
            stats.num_deriv_eval = iter;
        end

    end

else
    error('Unknown mode');
end