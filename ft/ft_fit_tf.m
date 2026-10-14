function [b, a, err, output] = ft_fit_tf(X, w, N_b, N_a, mode, options)
% Least squares fit of transfer function numerator / denominator coefficients
% to frequency response matrix X:

% min_{b, a} \sum_n wt(w_n) * | X_n - H_b(w_n) / H_a(w_n) |^2
%          = \sum_n wt(w_n) * | H_a(w_n) * X_n - H_b(w_n) |^2 / |H_a(w_n)|^2

% Iterative re-weighted least squares (Sanathanan-Koerner iterations):
% (b^{k}, a^{k})      = argmin_{b, a} \sum_n wt(w_n) / max(wt^{k}(w_n), epsilon) 
%                     * | H_a(w_n) * X_n - H_b(w_n) |^2

% b:        Numerator coefficients
% a:        Denominator coefficients
% wt:       Frequency weighting
% alpha:    Step-size between [0, 1]
% epsilon:  Regularization parameter

% Step-size weight controller: Update rule options.alpha between [0, 1] for wt_func
% wt^{k}(w_n)        = wt_func( H_a^{k-1}(w_n), wt^{k-1}(w_n), alpha )^2

% Transfer function types are specified in the mode input argument: 
% * Default is unconstrained coefficients b, a
% * If system is underdetermined, find the minimum norm ||[b, a]|| solution
% * minphase replaces target X with minimum phase reconstruction of real cepstrum
%   computed from interpolated |X|. Ignore stabilization constraints during iterations 
%   until final minimum phase allpass decomposition. 
% * allpass treats exp(1i * angle(x)) as target, constrains b = fliplr(a)

% Options
% * Weight function specification: See options.wt_mode
% * Stable denominator coefficients a: See options.stabilize_mode

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X:         [N x M] N frequency responses, M functions
%w:         [N x 1] Angular frequency (radians / sample) [0 to 2 * pi]
%N_b:       Number of numerator coefficients        \sum_{n=0)^N_b b_n * exp(-1i * w * n)
%N_a:       Number of denominator coefficients  1 + \sum_{n=1)^N_a a_n * exp(-1i * w * n)
%mode:      String, transfer function constraints {'unconstr', 'minphase', 'allpass'}
%                   'unconstr':     Unconstrained
%                   'minphase':     Poles-zeros within complex unit circle
%                   'allpass':      Conjugate symmetric poles-zeros

% options:              struct

% options.max_iter:     Maximum number of iterations
% options.wt_mode:      String, weight update mode { 'AM', 'GM', 'HM', 'BM'}
%                           'AM'    Arithmetic mean: wt_func(x,y) =  x * alpha + y * (1 - alpha) 
%                           'GM'    Geometric mean:  wt_func(x,y) =  x^alpha * y^(1 - alpha)
%                            HM'    Harmonic mean:   wt_func(x,y) =  1 / (alpha / x + (1 - alpha) / y)
%                           'BM'    Barycentric mean (inverse CDF interpolation)

% options.wt:           [N x 1]: Frequency weighting, [] for unity
% options.alpha:        Scalar, weight update rate between [0, 1]
% options.epsilon:      Scalar, minimum denominator weighting

% options.stabilize_mode:    String, stabilize poles  {'none', 'reflect', 'realz'}
%                               'none':     Allow unstable poles
%                               'reflect':  Reflect poles into complex unit circle per iteration
%                               'realz':    Enforce linear constrained real(H_b(w)) > 0 for grid of w in [0, pi]
%                                           Reference: A. Chottera and G. Jullien, "A linear programming approach to recursive digital filter design with linear phase," in IEEE Transactions on Circuits and Systems, vol. 29, no. 3, pp. 139-149, March 1982, doi: 10.1109/TCS.1982.1085123
% options.stabilize_M_w:     Number of uniform angular frequencies betwen [0, pi] for stabilization in 'stabilize_mode' = 'realz'
% options.stabilize_tol0:    Margin around unit diskin 'stabilize_mode' = 'realz

% options.enable_disp:  Logical, if true, plot fit
% options.Fs:           Sample rate for display
% options.disp_dB_lim:  [1 x 2] dB [min, max]
% options.disp_deg_lim: [1 x 2] degree [min, max]
% options.disp_output:  Logical, if true, include outputs in title display

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%b:         [1 x N_b]     Numerator coefficients
%a:         [1 x N_a]     Denominator coefficients
%err:       Scalar, wt' * abs((Fb * b(:) ./ (Fa * a(:))) - X).^2;

%output:    struct
%output.err_iters:  [options.max_iter x 1] Errors per iteration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Fit to randomized transfer function

% rng(51213);
% Fs = 48000;
% N_b_ref = 10;
% N_a_ref = 12;
% b_ref = randn(1, N_b_ref);
% a_ref = randn(1, N_a_ref);
% a_ref = polystab(a_ref);
% 
% freq = logspace(log10(20), log10(Fs/2), 512)';
% w = freq / Fs * (2 * pi);
% X = freqz(b_ref, a_ref, w);
% 
% N_b = 8;
% N_a = 10;
% [b_AM, a_AM] = ft_fit_tf(X, w, N_b, N_a, 'unconstr', 'wt_mode', 'AM', 'enable_disp', true);
% [b_GM, a_GM] = ft_fit_tf(X, w, N_b, N_a, 'unconstr', 'wt_mode', 'GM', 'enable_disp', true);
% [b_HM, a_HM] = ft_fit_tf(X, w, N_b, N_a, 'unconstr', 'wt_mode', 'HM', 'enable_disp', true);
% [b_BM, a_BM] = ft_fit_tf(X, w, N_b, N_a, 'unconstr', 'wt_mode', 'BM', 'enable_disp', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Fit to randomized minphase transfer function

% rng(5113);
% Fs = 48000;
% N_b_ref = 10;
% N_a_ref = 12;
% b_ref = randn(1, N_b_ref);
% a_ref = randn(1, N_a_ref);
% [b_ref, a_ref] = ft_tf_minphase(b_ref, a_ref);
% 
% freq = logspace(log10(20), log10(Fs/2), 512)';
% w = freq / Fs * (2 * pi);
% X = freqz(b_ref, a_ref, w);
% 
% N_b = 6;
% N_a = 8;
% [b_AM, a_AM] = ft_fit_tf(X, w, N_b, N_a, 'minphase', 'wt_mode', 'AM', 'enable_disp', true);
% [b_GM, a_GM] = ft_fit_tf(X, w, N_b, N_a, 'minphase', 'wt_mode', 'GM', 'enable_disp', true);
% [b_HM, a_HM] = ft_fit_tf(X, w, N_b, N_a, 'minphase', 'wt_mode', 'HM', 'enable_disp', true);
% [b_BM, a_BM] = ft_fit_tf(X, w, N_b, N_a, 'minphase', 'wt_mode', 'BM', 'enable_disp', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Fit to randomized allpass transfer function

% rng(51213);
% Fs = 48000;
% N_a_ref = 12;
% a_ref = randn(1, N_a_ref);
% a_ref = polystab(a_ref);
% b_ref = fliplr(a_ref);
% 
% freq = logspace(log10(20), log10(Fs/2), 512)';
% w = freq / Fs * (2 * pi);
% X = freqz(b_ref, a_ref, w);

% N_a = 10;
% [b_AM, a_AM] = ft_fit_tf(X, w, N_a, N_a, 'allpass', 'wt_mode', 'AM', 'enable_disp', true);
% [b_GM, a_GM] = ft_fit_tf(X, w, N_a, N_a, 'allpass', 'wt_mode', 'GM', 'enable_disp', true);
% [b_HM, a_HM] = ft_fit_tf(X, w, N_a, N_a, 'allpass', 'wt_mode', 'HM', 'enable_disp', true);
% [b_BM, a_BM] = ft_fit_tf(X, w, N_a, N_a, 'allpass', 'wt_mode', 'BM', 'enable_disp', true);

arguments
    X (:,1) double = 1;
    w (:,1) double = 0;

    N_b (1,1) double {mustBePositive, mustBeInteger} = 1;
    N_a (1,1) double {mustBePositive, mustBeInteger} = 1;

    mode (1,:) char {mustBeMember(mode, {'unconstr', 'minphase', 'allpass'} )} = 'unconstr';

    options.max_iter (1,1) double {mustBePositive, mustBeInteger} = 40;

    options.wt_mode  (1,:) char {mustBeMember(options.wt_mode, { 'AM', 'GM', 'HM', 'BM'})} = 'AM';
    options.wt (:,1) double {mustBePositive} = [];
    options.alpha (1,1) double {mustBeNonnegative} = 0.5;
    options.epsilon (1,1) double {mustBeNonnegative} = 1e-8;
    options.stabilize_mode (1,:) char {mustBeMember(options.stabilize_mode, {'none', 'reflect', 'realz'})} = 'none';
    options.stabilize_M_w (1,1) double {mustBePositive, mustBeInteger} = 128;
    options.stabilize_tol0 (1,1) double {mustBeNonnegative} = 1e-8; 

    options.enable_disp (1,1) logical = false;
    options.Fs (1,1) double {mustBePositive} = 48000;
    options.disp_dB_lim (1,2) double = [-inf, inf];
    options.disp_deg_lim (1,2) double = [-180, 180];
    options.disp_output (1,1) logical = true;
end

N = numel(X);
assert(N == numel(w), 'X, w size mismatch');

if isempty(options.wt)
    wt = ones(N, 1);
else
    wt = options.wt;
end
assert(N == numel(wt), 'X, options.wt size mismatch');
wt_sqrt = sqrt(wt);
epsilon_sqrt = sqrt(options.epsilon);

assert(options.alpha <= 1, 'options.alpha must between [0, 1]');

% Denominator weight function
if strcmp(options.wt_mode, 'AM')

    wt_func = @(x, y)  (x .* options.alpha) + (y .* (1 - options.alpha));

elseif strcmp(options.wt_mode, 'GM')

    wt_func = @(x, y)  (x.^options.alpha) .* (y.^(1 - options.alpha));

elseif strcmp(options.wt_mode, 'HM')

    wt_func = @(x, y) 1 ./ ( (options.alpha ./ x) + (1 - options.alpha) ./ y);

elseif strcmp(options.wt_mode, 'BM')
    
    wt_func = @(x, y) cdf_bm(x, y, options.alpha);

else
    error('Unsupported options.wt_mode');
end

% Pre-compute DFT matrices
Fb = exp(-1i * w * (0:(N_b-1)) );  
Fa = exp(-1i * w * (0:(N_a-1)) );

% Initial values and states
H_a_mag_prev = ones(N, 1);
wt_func_prev = ones(N, 1); 
b = [1, zeros(1, N_b-1)];
a = [1, zeros(1, N_a-1)];
err = wt' * abs((Fb * b(:) ./ (Fa * a(:))) - X).^2;
err_iters = zeros(options.max_iter, 1);

if strcmp(mode, 'unconstr') 

    % Pre-compute stabilization constraints
    if strcmp(options.stabilize_mode, 'realz')
    
        Fa_stab = real(exp(-1i * linspace(0, pi, options.stabilize_M_w)' * (0:(N_a-1)) ));
        A_stab_ineq = -[Fa_stab(:, 2:end), zeros(options.stabilize_M_w, N_b)];
        b_stab_ineq = (1 - options.stabilize_tol0) * Fa_stab(:, 1);
        
        options_quadprog = optimoptions('quadprog', 'Display', 'none');
        
    end

    for k = 1:options.max_iter

        wt_func_curr    = wt_func(H_a_mag_prev, wt_func_prev);
        wt_func_prev    = wt_func_curr;        
        wt_sqrt_k       = wt_sqrt ./ max(wt_func_curr, epsilon_sqrt);
        
        %figure; semilogx(w/(2*pi)*options.Fs, wt_func_curr, w/(2*pi)*options.Fs, H_a_mag_prev, 'linewidth', 2); grid on; axis tight; legend('curr', 'H', 'location', 'best');

        F = bsxfun(@times, [bsxfun(@times, Fa(:, 2:end), X), -Fb], wt_sqrt_k);
        Y = -(wt_sqrt_k .* X);

        if strcmp(options.stabilize_mode, 'none')

            x = lsqminnorm([real(F); imag(F)], [real(Y); imag(Y)]);
            b_k = x(N_a:end).';
            a_k = [1; x(1:(N_a - 1))].';

        elseif strcmp(options.stabilize_mode, 'reflect')

            x = lsqminnorm([real(F); imag(F)], [real(Y); imag(Y)]);
            b_k = x(N_a:end).';
            a_k = [1; x(1:(N_a - 1))].';
            a_k = polystab(a_k); 

        elseif strcmp(options.stabilize_mode, 'realz')

            F_re_im = [real(F); imag(F)];
            Y_re_im = [real(Y); imag(Y)];

            if size(F_re_im, 1) >= size(F_re_im, 2) % Exact or overdetermined 
                
                [x, fval, exitflag] = quadprog(F_re_im' * F_re_im, -F_re_im' * Y_re_im, A_stab_ineq, b_stab_ineq, ...
                    [], [], [], [], [], options_quadprog);

            else % Underdetermined

                % Try to find minimum norm solution
                 [x, fval, exitflag] = quadprog(eye(N_b + N_a - 1), zeros(N_b + N_a - 1, 1), A_stab_ineq, b_stab_ineq, ...
                     [real(F); imag(F)], [real(Y); imag(Y)], [], [], [], options_quadprog);
                 if exitflag ~= 1 % Revert to exact or overdetermined case
                    [x, fval, exitflag] = quadprog(F_re_im' * F_re_im, -F_re_im' * Y_re_im, A_stab_ineq, b_stab_ineq, ...
                        [], [], [], [], [], options_quadprog);
                 end
            end

            b_k = x(N_a:end).';
            a_k = [1; x(1:(N_a - 1))].';

            assert(all(abs(roots(a_k)) <= 1), 'unstable poles');
            
        else
            error('Unsupported options.stabilize_mode');
        end

        H_a_k = Fa * a_k(:);
        H_a_mag_prev = abs(H_a_k);

        % Compute error
        H_b_k = Fb * b_k(:);
        H_k = H_b_k ./ H_a_k;
        err_iters(k) = wt' * abs(H_k - X).^2;
        if err_iters(k) < err
            err = err_iters(k);
            b = b_k;
            a = a_k;
        end
    
    end

elseif strcmp(mode, 'minphase') 

    [~, h_mp_rceps] = ft_fit_minphase(abs(X), w, 1);
    X = freqz(h_mp_rceps, 1, w);

    for k = 1:options.max_iter

        wt_func_curr    = wt_func(H_a_mag_prev, wt_func_prev);
        wt_func_prev    = wt_func_curr;        
        wt_sqrt_k       = wt_sqrt ./ max(wt_func_curr, epsilon_sqrt);
        
        %figure; semilogx(w/(2*pi)*options.Fs, wt_func_curr, w/(2*pi)*options.Fs, H_a_mag_prev, 'linewidth', 2); grid on; axis tight; legend('curr', 'H', 'location', 'best');

        F = bsxfun(@times, [bsxfun(@times, Fa(:, 2:end), X), -Fb], wt_sqrt_k);
        Y = -(wt_sqrt_k .* X);
        x = lsqminnorm([real(F); imag(F)], [real(Y); imag(Y)]);
        
        b_k = x(N_a:end).';
        a_k = [1; x(1:(N_a - 1))].';
                       
        H_a_k = Fa * a_k(:);
        H_a_mag_prev = abs(H_a_k);
     
        % Compute error
        H_b_k = Fb * b_k(:);
        H_k = H_b_k ./ H_a_k;
        err_iters(k) = wt' * abs(H_k - X).^2;
        if err_iters(k) < err
            err = err_iters(k);
            b = b_k;
            a = a_k;
        end
    
    end

    if ~isminphase(b, a)
        [b, a] = ft_tf_minphase(b, a);
        err = wt' * abs(freqz(b, a, w) - X).^2;
    end

elseif strcmp(mode, 'allpass')

    assert(N_b == N_a, 'N_b must equal N_a');
    X = X ./ abs(X);

    % Pre-compute stabilization constraints
    if strcmp(options.stabilize_mode, 'realz')
    
        Fa_stab = real(exp(-1i * linspace(0, pi, options.stabilize_M_w)' * (0:(N_a-1)) ));
        A_stab_ineq = -Fa_stab(:, 2:end);
        b_stab_ineq = (1 - options.stabilize_tol0) * Fa_stab(:, 1);

        options_quadprog = optimoptions('quadprog', 'Display', 'none');
        
    end

    for k = 1:options.max_iter

        wt_func_curr    = wt_func(H_a_mag_prev, wt_func_prev);
        wt_func_prev    = wt_func_curr;        
        wt_sqrt_k       = wt_sqrt ./ max(wt_func_curr, epsilon_sqrt);
        
        %figure; semilogx(w/(2*pi)*options.Fs, wt_func_curr, w/(2*pi)*options.Fs, H_a_mag_prev, 'linewidth', 2); grid on; axis tight; legend('curr', 'H', 'location', 'best');

        F = bsxfun(@times, bsxfun(@times, Fa(:, 2:end), X) - fliplr(Fa(:, 1:end-1)), wt_sqrt_k);
        Y = wt_sqrt_k .* (Fa(:, end) - X);
        if strcmp(options.stabilize_mode, 'none')

            x = lsqminnorm([real(F); imag(F)], [real(Y); imag(Y)]);
            a_k = [1; x(1:(N_a - 1))].';
            b_k = fliplr(a_k);

        elseif strcmp(options.stabilize_mode, 'reflect')

            x = lsqminnorm([real(F); imag(F)], [real(Y); imag(Y)]);
            a_k = [1; x(1:(N_a - 1))].';
            a_k = polystab(a_k);
            b_k = fliplr(a_k);

        elseif strcmp(options.stabilize_mode, 'realz')

            F_re_im = [real(F); imag(F)];
            Y_re_im = [real(Y); imag(Y)];

            if size(F_re_im, 1) >= size(F_re_im, 2) % Exact or overdetermined 
                
                [x, fval, exitflag] = quadprog(F_re_im' * F_re_im, -F_re_im' * Y_re_im, A_stab_ineq, b_stab_ineq, ...
                    [], [], [], [], [], options_quadprog);

            else % Underdetermined

                % Try to find minimum norm solution
                 [x, fval, exitflag] = quadprog(eye(N_a - 1), zeros(N_a - 1, 1), A_stab_ineq, b_stab_ineq, ...
                     [real(F); imag(F)], [real(Y); imag(Y)], [], [], [], options_quadprog);
                 if exitflag ~= 1 % Revert to exact or overdetermined case
                    [x, fval, exitflag] = quadprog(F_re_im' * F_re_im, -F_re_im' * Y_re_im, A_stab_ineq, b_stab_ineq, ...
                        [], [], [], [], [], options_quadprog);
                 end
            end

            a_k = [1; x(1:(N_a - 1))].';
            b_k = fliplr(a_k);

            assert(all(abs(roots(a_k)) <= 1), 'unstable poles');
            
        else
            error('Unsupported options.stabilize_mode');            
        end
                               
        H_a_k = Fa * a_k(:);
        H_a_mag_prev = abs(H_a_k);
     
        % Compute error
        H_b_k = Fb * b_k(:);
        H_k = H_b_k ./ H_a_k;
        err_iters(k) = wt' * abs(H_k - X).^2;
        if err_iters(k) < err
            err = err_iters(k);
            b = b_k;
            a = a_k;
        end
    
    end

else 

    error('Unsupported mode');

end

% Additional outputs
output = [];
output.err_iters = err_iters;
output.stable = isstable(b, a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.enable_disp && coder.target("MATLAB")

    fontsize = 16;
    linewidth = 2;

    hz_disp = logspace(log10(20), log10(options.Fs/2), 512);
    H_disp = freqz(b, a, hz_disp, options.Fs);

    hz_in = w / (2*pi) * options.Fs;
    H_in = freqz(b, a, hz_in, options.Fs);

    h_fig = figure;
    h_fig.Position = [100, 100, 800, 600];
    h_t = tiledlayout(2, 1);

    title_str = ['Transfer Function Fit: Mode = ', mode, 'Numer = ', num2str(N_b), ', Denom = ', num2str(N_a), ', Update = ', options.wt_mode];
    if options.disp_output
         title_str = [title_str, '\newlineError = ', num2str(err), ', Stable = ', num2str(output.stable)];
    end
    title(h_t, title_str, 'fontsize', fontsize + 2);

    % Magnitude
    nexttile;
    semilogx(hz_in, mag2db(abs(X)), 'ro-', 'markersize', 4, 'linewidth', linewidth); hold on;
    semilogx(hz_disp, mag2db(abs(H_disp)), 'b-', hz_in, mag2db(abs(H_in)) - mag2db(abs(X)), 'g-', 'linewidth', 1.5 * linewidth); hold on;

    yline(0, 'k:', 'linewidth', 0.8 * linewidth);
    h_lg =  legend('Target Response', 'Fitted Response', 'Residual', ...
        'location', 'best', 'Orientation', 'horizontal');  set(h_lg, 'fontsize', fontsize - 1);
    xlabel('Frequency (Hz)', 'fontsize', fontsize); 
    ylabel('Magnitude (dB)', 'fontsize', fontsize);
    title('Magnitude Response', 'fontsize', fontsize + 1);
    grid on; axis tight;
    set(gca, 'fontsize', fontsize - 1);
    ylim(options.disp_dB_lim);
    xlim([min(hz_disp), max(hz_disp)]);

    % Phase
    nexttile;
    semilogx(hz_in, rad2deg(angle(X)), 'ro-', 'markersize', 4, 'linewidth', linewidth); hold on;
    semilogx(hz_disp, rad2deg(angle(H_disp)), 'b-', hz_in, rad2deg(angle(H_in ./ X)), 'g-', 'linewidth', 1.5 * linewidth); hold on;

    yline(0, 'k:', 'linewidth', 0.8 * linewidth);
    h_lg =  legend('Target Response', 'Fitted Response', 'Residual', ...
        'location', 'best', 'Orientation', 'horizontal');  set(h_lg, 'fontsize', fontsize - 1);
    xlabel('Frequency (Hz)', 'fontsize', fontsize); 
    ylabel('Phase (Degrees)', 'fontsize', fontsize);
    title('Phase Response', 'fontsize', fontsize + 1);
    grid on; axis tight;
    set(gca, 'fontsize', fontsize - 1);
    ylim(options.disp_deg_lim);
    xlim([min(hz_disp), max(hz_disp)]);

    err_iters
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Barycentric weighted mean from 1D cumulative distributions of f, g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%f:     [N x 1]
%g:     [N x 1]
%alpha: Scalar between [0, 1]

%options:   struct
%options.M: Number of CDF interpolation points for u

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h:     [N x 1]

function h = cdf_bm(f, g, alpha, options)

arguments
    f (:,1) double {mustBePositive} = 1;
    g (:,1) double {mustBePositive} = 1;
    alpha (1,1) double {mustBeNonnegative} = 1;

    options.M (1,1) double {mustBePositive} = 4096;
end

N = numel(f);
assert(N == numel(g), 'f, g size mismatch');
assert(alpha <= 1, 'alpha must be between [0, 1]');

if alpha == 0

    h = y;

elseif alpha == 1

    h = x;

else

    sum_f = sum(f);
    F = [0; cumsum(f)] / sum_f;
    F(end) = 1;

    sum_g = sum(g);
    G = [0; cumsum(g)] / sum_g;
    G(end) = 1;

    w = linspace(0, pi, N + 1)';
    u = linspace(0, 1, options.M)'; % Sampling points along inverse CDF

    F_inv = interp1(F, w, u); 
    G_inv = interp1(G, w, u);
    H_inv = alpha * F_inv + (1 - alpha) * G_inv;

    H = interp1(H_inv, u, w);
    h = diff(H) * (sum_f * alpha + sum_g * (1 - alpha) );
 
    ;
    %figure;  yyaxis left; plot(w(1:end-1), f, 'r-', w(1:end-1), g, 'b-', w(1:end-1), h, 'k-', 'linewidth', 2);
    % yyaxis right; plot(w, F, 'r--', w, G, 'b--', w, H, 'k--', 'linewidth', 2); grid on; axis tight; legend('f', 'g', 'h', 'F', 'G', 'H', 'location', 'best');

end
  
