function [g_exp, G_exp] = ft_exp_conv_g(g, n, max_N)
%Compute convolution of g with itself n times with memoization
%if max_N = 0, clear persistent variables

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%g:         [N x 1] Exponentiating IR
%n:         Number of self-convolutions
%max_N:     Number of samples in output spectrum, clear persistent if 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%g_exp:     [min(n * (N - 1) + 1, max_N) x 1]   Self-convolution
%G_exp:     [max_N x 1]                         Spectrum of self-convolutions


persistent g_n      %Cell-array, g_n{n} = g convolved with itself n times
persistent G_n      %Cell-array, G_n{n} = spectrum of g convolved with itself n times

if max_N == 0
    clear g_n;
    clear G_n;
    return;
end

N = numel(g);

%Expand cell array
if n >= numel(g_n)
    g_n{n} = [];
    G_n{n} = [];
end

%Base-case
if n == 1
    g_n{1} = g;
    G_n{1} = fft(g, max_N);
end

if isempty(g_n{n})
    n_half = floor(n / 2);
    [g_low,     G_low]     = ft_exp_conv_g(g, n_half, max_N);
    [g_high,    G_high]    = ft_exp_conv_g(g, n - n_half, max_N);
 
    G_n{n} = G_low .* G_high;
    g_n{n} = ifft(G_n{n});
    g_n{n} = g_n{n}(1:min(n * (N - 1) + 1, numel(g_n{n}) ));

end
g_exp = g_n{n};
G_exp = G_n{n};