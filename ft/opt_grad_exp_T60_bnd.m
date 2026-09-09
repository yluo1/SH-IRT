function [m_max, y_max, y_max_apx] = opt_grad_exp_T60_bnd(T60, Fs, enable_disp)
%Find argmax_m y(m)  = m * |G|^(m-1) for G = 10^(-3 / (Fs * T60) )
%for the gradient of the time-varying filter's magnitude under exponentiation
%w.r.t. T60(z, theta, phi) and sample rate Fs

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%T60:           60 dB attenuation time (seconds)
%Fs:            Sample rate
%enable_disp:   Logical, if true, display y(m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%m_max:         Solution m
%y_max:         Maximal value
%y_max_apx:     Linear approximation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample Usage:

%[m_max, y_max, y_max_apx] = opt_grad_exp_T60_bnd(1.5, 48000, true);

arguments
    T60 (1,1) double {mustBePositive} = 1;
    Fs (1,1) double {mustBePositive} = 48000;
    enable_disp (1,1) logical = false;
end

syms m_sym Fs_sym T_sym
y_sym = m_sym .* 10.^(-3 * (m_sym - 1) ./ (Fs_sym * T_sym))

dy_dm = simplify(diff(y_sym, m_sym), 100)

m_sol = solve(dy_dm == 0)

y_max = simplify(subs(y_sym, m_sym, m_sol), 100)


m_max = double(subs(subs(m_sol, Fs_sym, Fs), T_sym, T60));
y_max = double(subs(subs(subs(y_max, Fs_sym, Fs), T_sym, T60), m_sym, m_max));

y_max_num = 10^(3/(Fs * T60)) * Fs * T60 * exp(-1) / (3*log(10));
y_max_apx = Fs * T60 / (3 * log(10) * exp(1));

% Equivalences
% norm(y_max - y_max_num)
% norm(y_max - (m_max *  10^( (m_max^(-1) - 1) / log(10) ) ) )
% norm(y_max - (m_max *  exp(1)^(m_max^(-1) - 1)  ))

%Plotting
if enable_disp

    m = 1:100000;
    y = m .* 10.^(-3 * (m - 1) ./ (Fs * T60));

    fontsize = 14;
    figure;
    semilogx(m, y, m_max, y_max, 'r*', m_max, y_max_apx, 'bo', 'linewidth', 1.5);
    xlabel('Sample m', 'fontsize', fontsize);
    ylabel('y(m)', 'fontsize', fontsize);
    set(gca, 'fontsize', fontsize - 1);
    h_lg = legend('y', 'max', 'linear approx max', 'location', 'best');
    title('Gradient Filter Response Magnitude: $y(m) = m |G|^{m-1}$', 'fontsize', fontsize + 3, 'interpreter', 'latex');
    set(h_lg, 'fontsize', fontsize - 1);
    
    grid on; axis tight;
end