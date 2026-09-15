function [f, h_fig] = ft_rec_conv_opt(h, g, options)
%Recursive convolution of FIR h by exponentiating FIR g

%Author: Yuancheng Luo, 2026

%Paper Reference:
%Yuancheng Luo, "Fast Time-Varying Exponentiated Convolution Methods for Generative Direction Dependent Reverberation",
%Proceedings of the 161th Audio Engineering Society Convention.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%h:                 [1 x M] Fixed IR
%g:                 [1 x N] Exponentiating IR

%options:           struct

%options.mode:      String, compute method 'fourblock', 'twoblock'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%f:                 [1 x M * N]     Filtered IR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Recursive convolution with two and three tap filters

% rng(4123);
% Fs = 48000;
% T = 1;
% M = ceil(T * Fs);
% h = randn(1, M);
% g = ft_two_tap_FIR(0, -1);
% f_two_tap     = ft_rec_conv_opt(h, g, 'enable_disp', true);
% f_three_tap   = ft_rec_conv_opt(h, [g, -0.0001], 'enable_disp', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare two block and four block subdivision modes

% rng(4123);
% M = 3777;
% h = randn(1, M);
% g = [ft_two_tap_FIR(0, -1), -0.01];
% tic; f_fourblock = ft_rec_conv_opt(h, g, 'mode', 'fourblock'); toc_fourblock = toc
% tic; f_twoblock  = ft_rec_conv_opt(h, g, 'mode', 'twoblock');  toc_twoblock = toc
% err = norm(f_fourblock - f_twoblock)

arguments
    h (1,:) double = [1 0];
    g (1,:) double = 1;
    
    options.mode {mustBeMember(options.mode, {'fourblock', 'twoblock'})} = 'twoblock';

    %Display options
    options.enable_disp (1,1) logical = false;
    options.disp_mode {mustBeMember(options.disp_mode, {'RIR', 'spec'}) } = 'spec'
    options.disp_trunc (1,1) logical = true;
    options.disp_position (1,4) double = [100, 100, 800, 400];

    options.N_FFT (1,1) double {mustBePositive, mustBeInteger} = 4096;
    options.Fs (1,1) double {mustBePositive, mustBeInteger} = 48000;

    options.dB_lim (1,2) double = [-80, inf];
    options.time_lim (1,2) double = [-inf, inf];
    options.freq_lim (1,2) double = [-inf, inf];

    options.fontsize = 18;
end

M = numel(h);

if strcmp(options.mode, 'fourblock')
    f = exp_conv_recur_fourblock(h(:).', g(:).');
else    %'twoblock'
    f = exp_conv_recur_twoblock(M, h(:).', g(:).');
end

%Plotting
h_fig = [];
if options.enable_disp && coder.target('MATLAB')
    fontsize = options.fontsize;
    
    if options.disp_trunc
        f_disp = f(1:M);
    else
        f_disp = f;
    end

    if strcmp(options.disp_mode, 'spec')
        
        h_fig = figure;
        h_fig.Position = options.disp_position;
        tiledlayout(1, 2);
        nexttile;
        spectrogram(h, options.N_FFT, ceil(options.N_FFT * 0.9), options.N_FFT, options.Fs, "power", "yaxis");
        axis tight; set(gca, 'YScale', 'log'); set(gca, 'fontsize', fontsize - 1);
        cb = colorbar; ylabel(cb, 'Power (dB)','FontSize', fontsize); title('Input Impulse Response', 'fontsize', fontsize + 1);
        colormap(hot);
        clim(options.dB_lim);
        xlim(options.time_lim);
        ylim(options.freq_lim/1000);

        nexttile;
        spectrogram(f_disp, options.N_FFT, ceil(options.N_FFT * 0.9), options.N_FFT, options.Fs, "power", "yaxis"); 
        axis tight; set(gca, 'YScale', 'log');  set(gca, 'fontsize', fontsize - 1);
        cb = colorbar; ylabel(cb, 'Power (dB)','FontSize', fontsize); title('Output Impulse Response', 'fontsize', fontsize + 1);
        colormap(hot);
        clim(options.dB_lim);
        xlim(options.time_lim);
        ylim(options.freq_lim/1000);

    
    elseif strcmp(options.disp_mode, 'RIR')

        options_plot_RIR = plot_RIR_opts('Fs', options.Fs);
        h_fig = plot_RIR(f_disp, options_plot_RIR);

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Recursive formulation for all sizes h, g, half-tail
%m:
%h: [1 x M]
%g: [1 x N]

%Output
%f:     [1 x M + m * (N - 1)]

function f = exp_conv_recur_twoblock(m, h, g)

M = numel(h);
N = numel(g);

if m <= 32
    %Problem size small, do direct computation
    f = exp_conv_m(m, h, g);
else
    %Problem size large
    m_half = floor(m / 2);

    %High block
    m_conv = m_half * (N - 1) + (M - m_half);
    G = fft(g, m_conv).^m_half;
    H = fft(h((m_half + 1):M), m_conv);
    
    f_hi = ifft(H .* G);

    %Low block
    f_lo = exp_conv_recur_twoblock(m_half, h(1:m_half), g); %[1 x m_half * N]
    
    %Combine blocks
    m_half_bar = M +  m_half * (N - 1);
    if coder.target('MATLAB')
        f_m_half = zeros(1, m_half_bar);
    else
        f_m_half = complex(zeros(1, m_half_bar));
    end
    f_m_half(1:(N * m_half))    = f_m_half(1:(N * m_half)) + f_lo;
    f_m_half((m_half+1):m_half_bar)  = f_m_half((m_half+1):m_half_bar) + f_hi;

    %Recur on truncated f_m_half
    f = [f_m_half(1:m_half), exp_conv_recur_twoblock(m - m_half, f_m_half((m_half+1):end), g)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Direct version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%m
%h:         [1 x M]
%g:         [1 x N]

%Output
%f:         [1 x M + m * (N - 1)]

function f = exp_conv_m(m, h, g)
M = numel(h);
N = numel(g);

size_f = M + m * (N - 1);

if coder.target('MATLAB')
    f = zeros(1, size_f);
else
    f = complex(zeros(1, size_f));
end
f(1:M) = h;
for i = 1:m
    idx_in      = i : (M + (i - 1) * (N - 1) );
    idx_out     = i : (M + i * (N - 1) );
    f(idx_out)  = conv(f(idx_in), g);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Recursive formulation for all sizes h, g, 4 blocks
%h:         [1 x M]
%g:         [1 x N]

function y = exp_conv_recur_fourblock(h, g)

M = numel(h);
N = numel(g);

y = zeros(1, M * N);
if M <= 32
    %Problem size small, do direct computation
    y = exp_conv(h, g);
else
    %Problem size large
    M_half = floor(M/2);
    
    %Left block
    y_low = exp_conv_recur_fourblock(h(1 : M_half), g);
%   y_low = exp_conv(h(1 : M_half), g);
    y(1:M_half) = y_low(1:M_half);
    
    %Mid Block
    M_mid = M - (M_half + 1) +1; %Size of the mid block
    M_pow = M_half;
    M_mid_padded = M_mid + M_pow * (N-1);  
    G = fft(g, M_mid_padded);
    H = fft( h((M_half+1) : M), M_mid_padded);
    y_mid = ifft((G .^ M_pow) .* H);
    
	%Top block
    buf = y_mid(1:M_mid);   
    if N*M_half < M
        buf(1:M_half) = buf(1:M_half) + y_low((M_half + 1):(N * M_half));
    else
        buf = buf + y_low((M_half + 1):M);
    end    

    y_top = exp_conv_recur_fourblock(buf, g);
%   y_top = exp_conv(buf, g);
    y((M_half + 1) : (M_half + N * M_mid)) = y_top;
    
    %Right block
    M_pad_2 = M*(N-1);
    buf = y_mid((M_mid+1) : M_mid_padded);
    buf(1:(N * M_half - M)) = buf(1:(N * M_half - M)) + y_low((M + 1) : N * M_half);
    R = fft(buf, M_pad_2);
    H2 = fft(g, M_pad_2);
    y_right = ifft((H2 .^ (M_mid) ) .* R);
	y((M+1) : (M + M_pad_2 ) ) =  y((M + 1) : (M + M_pad_2 ) ) + y_right;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Direct version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%h:         [1 x M]
%g:         [1 x N]

function y = exp_conv(h, g)
M = numel(h);
N = numel(g);

size_y = M * N;

if N == 2
    %Unrolled 2-tap reverse order variant
    y = zeros(1, size_y);
    y(1:M) = h;
    for i = 1:M
        for j = size_y : -1 : (i+1)
            y(j) = y(j) * g(1) + y(j-1) * g(2);
        end
        y(i) = y(i) * g(1);
    end
else
    y = zeros(1, size_y);
    y(1:M) = h;
    for i = 1:M
        idx_in      = i : (M + (i - 1) * (N - 1) );
        idx_out     = i : (M + i * (N - 1) );
        y(idx_out)  = conv(y(idx_in), g);
    end

end




