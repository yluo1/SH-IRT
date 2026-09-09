function [f, h_fig] = ft_exp_conv_opt(h, g, options)
%Convolution of FIR h by time-varying exponentiating FIR g

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%h:                 [1 x M] Fixed IR
%g:                 [1 x N] Exponentiating IR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%f:                 [1 x M * N]     Filtered IR
%h_fig:             Figure handle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:

% rng(4123);
% Fs = 48000;
% T = 1;
% M = ceil(T * Fs);
% h = randn(1, M);
% g = ft_two_tap_FIR(0, -1);
% f = ft_exp_conv_opt(h, g, 'enable_disp', true);
% f2 = ft_exp_conv_opt(h, [g, -0.0001], 'enable_disp', true);

arguments
    h (1,:) double = [1 0];
    g (1,:) double = 1;
    
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

ft_exp_conv_g([], [], 0); %Clear persistent variables

M = numel(h);
N = numel(g);

h = h(:).';
g = g(:).';

if M == 1
    f = h * g;
else
    f = exp_conv_recur(h, g, M * (N - 1) + 1 );
end
%f_ref = exp_conv(h, g);

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
        cb = colorbar;
        ylabel(cb, 'Power (dB)','FontSize', fontsize); title('Input Impulse Response', 'fontsize', fontsize + 1);
        colormap(hot);
        clim(options.dB_lim);
        xlim(options.time_lim);
        ylim(options.freq_lim/1000);

        nexttile;
        spectrogram(f_disp, options.N_FFT, ceil(options.N_FFT * 0.9), options.N_FFT, options.Fs, "power", "yaxis"); 
        axis tight; set(gca, 'YScale', 'log');  set(gca, 'fontsize', fontsize - 1);
        cb = colorbar;
        ylabel(cb, 'Power (dB)','FontSize', fontsize); title('Output Impulse Response', 'fontsize', fontsize + 1);
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
%Direct version. Use for small size M

%Input
%h:         [1 x M] Fixed IR
%g:         [1 x N] Exponentiating IR

%Output
%f:         [1 x M * N] IR 

function f = exp_conv(h, g)
M = numel(h);
N = numel(g);
f = zeros(1, N * M);
g_orig = g;
for m = 1:M
    idx = m:(N * m);
    f(idx) = f(idx) + h(m) * g;
    g = conv(g, g_orig);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Recurrent version. F(x, h, M) = F(x(1:M/2), h, M/2) + delay(M/2+1) * conv(F(x(M/2+1:M), h, M-N), H^(M/2) )

%Input
%h:         [1 x M] Fixed IR
%g:         [1 x N] Exponentiating IR 
%max_N:     Max FFT size

%Output
%f:         [1 x M*N]

function f = exp_conv_recur(h, g, max_N)

M = numel(h);
M_half = floor(M / 2);
N = numel(g);
MN = M * N;

if coder.target('MATLAB')
    f = zeros(1, MN);
else
    f = complex(zeros(1, MN));
end

if M_half <= 32
    f_low = exp_conv(h(1:M_half), g); 
else
    f_low = exp_conv_recur(h(1:M_half), g, max_N); 
end

if M - M_half <= 32
    f_high = exp_conv(h(M_half+1:M), g);
else
    f_high = exp_conv_recur(h(M_half+1:M), g, max_N);
end

%Compute in frequency domain with memoization
%g_M_half = ft_exp_conv_g(g, M_half, max_N);
%f_high = conv(f_high, g_M_half);

%Equivalent without memoization
%f_high = ifft(fft(f_high, MN - M_half) .* fft(g_M_half, MN - M_half));  
f_high = ifft(fft(f_high, MN - M_half) .* (fft(g, MN - M_half).^M_half));

%Sum halves
idx_low     = 1:(M_half * N);
f(idx_low)  = f(idx_low) + f_low;

idx_high    = (M_half + 1):MN;
f(idx_high) = f(idx_high) + f_high;
