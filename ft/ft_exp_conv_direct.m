function f = ft_exp_conv_direct(h, g, mode)
%Time-varying exponentiated convolutions. %Direct version, use for small size M.

%Author: Yuancheng Luo, 2026

%Paper Reference:
%Yuancheng Luo, "Fast Time-Varying Exponentiated Convolution Methods for Generative Direction Dependent Reverberation",
%Proceedings of the 161th Audio Engineering Society Convention.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%h:         [1 x M] Fixed IR
%g:         [1 x N] Exponentiating IR

%mode:      String, compute mode
%           'recur'     Recurence relation (slow)
%           'conv'      Convolution equivalent

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%f:         [1 x M * N] Full IR 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare conv with recur modes

% M = 20;
% h = randn(1, M);
% g = ft_two_tap_FIR(0, -1);
% tic; f_conv = ft_exp_conv_direct(h, g, 'conv'); toc_conv = toc
% tic; f_recur = ft_exp_conv_direct(h, g, 'recur'); toc_recur = toc

% err = norm(f_conv - f_recur)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare with optimized FFT variant

% rng(4123);
% M = 10000;
% h = randn(1, M);
% g = [ft_two_tap_FIR(0, -1), -0.001];
% tic; f_conv        = ft_exp_conv_direct(h, g, 'conv'); toc_conv = toc
% tic; [~, f_fft]    = ft_exp_conv_opt(h, g); toc_fft = toc
% err                = norm(f_conv(:) - f_fft(:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare relative spectrum

% rng(4123);
% M = 200;
% h = randn(1, M);
% h2 = randn(1, M);
% g = ft_two_tap_FIR(0, -1);

% f = ft_exp_conv_direct(h, g);
% f2 = ft_exp_conv_direct(h2, g);

% %Plotting
% N_freq = 256;
% Fs = 48000;
% freq = logspace(log10(20), log10(Fs / 2) );
% H     = freqz(h, 1, freq, Fs);
% H2    = freqz(h2, 1, freq, Fs);
% F     = freqz(f, 1, freq, Fs);
% F2    = freqz(f2, 1, freq, Fs);

% plot_cpx_resp([H(:) ./ H2(:), F(:) ./ F2(:)], freq, 'resp_names', {'H(\omega) / H2(\omega)', 'F(\omega) / F2(\omega)'})

arguments
    h (1,:) double = [1 0];
    g (1,:) double = 1;
    mode {mustBeMember(mode, {'conv', 'recur'}) } = 'conv';
end

M = numel(h);
N = numel(g);
f = zeros(1, N * M);

if strcmp(mode, 'recur')
    
    MN = M * N;
    for n = 1:MN
        for m = 1:MN

            idx = n - m + 1;
            if idx >= 1 && idx <= M
                f(n) = f(n) + h(idx) * g_recur(g, idx, m);
            end

        end
    end
    
else %'conv'

    g_orig = g;
    for m = 1:M
        idx = m:(N * m);
        f(idx) = f(idx) + h(m) * g;
        g = conv(g, g_orig);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g_m_n = g_recur(g, m, n)

N = numel(g);

if m == 1 %Base case

    if n >= 1 && n <= N
        g_m_n = g(n);
    else
        g_m_n = 0;
    end

elseif m > 1 %Higher-order conv

    g_m_n = 0;
    for k = 1:(m * N)
        idx = n - k + 1;
        if idx >= 1 && idx <= N
            g_m_n = g_m_n + g(idx) * g_recur(g, m - 1, k);
        end
    end

else %m < 1

    g_m_n = 0;

end