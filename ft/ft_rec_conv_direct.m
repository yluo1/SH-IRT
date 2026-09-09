function f = ft_rec_conv_direct(h, g, mode)
%Recursive convolution of h by g. Direct version, use for small size M

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%h:         [1 x M] Fixed IR
%g:         [1 x N] Exponentiating IR

%mode:      String, compute mode
%           'conv'
%           'recur'
%           'recur_swap'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%f:         [1 x M * N] Full IR 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare conv with recur modes

% rng(4123);
% M = 20;
% h = randn(1, M);
% g = ft_two_tap_FIR(0, -1);
% tic; f_conv = ft_rec_conv_direct(h, g, 'conv'); toc_conv = toc
% tic; f_recur = ft_rec_conv_direct(h, g, 'recur'); toc_recur = toc
% tic; f_recur_swap = ft_rec_conv_direct(h, g, 'recur_swap'); toc_recur_swap = toc

% err_recur         = norm(f_conv - f_recur)
% err_recur_swap    = norm(f_conv - f_recur_swap)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare with optimized FFT variant

% rng(4123);
% M = 10000;
% h = randn(1, M);
% g = [ft_two_tap_FIR(0, -1), -0.001];
% tic; f_conv        = ft_rec_conv_direct(h, g, 'conv'); toc_conv = toc
% tic; [f_fft]       = ft_rec_conv_opt(h, g); toc_fft = toc
% err                = norm(f_conv(:) - f_fft(:))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare relative H(w)/H2(w) versus F(w)/F2(w)

% rng(4123);
% M = 200;
% h = randn(1, M);
% h2 = randn(1, M);
% g = ft_two_tap_FIR(0, -1);

% f = ft_rec_conv_direct(h, g);
% f2 = ft_rec_conv_direct(h2, g);

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

    mode {mustBeMember(mode, {'conv', 'recur', 'recur_swap'}) } = 'conv'
end

M = numel(h);
N = numel(g);

MN = M * N;
f = zeros(1, MN);

if strcmp(mode, 'recur')
    
    for n = 1:MN
        f(n) = f_recur(h, g, M, n);
    end

elseif strcmp(mode, 'recur_swap')
    
    for n = 1:MN
        f(n) = f_recur_swap(h, g, M, n);
    end

else %'conv'

    f(1:M) = h;
    for m = 1:M
        idx_in      = m : (M + (m - 1) * (N - 1) );
        idx_out     = m : (M + m * (N - 1) );
        f(idx_out)  = conv(f(idx_in), g);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = f_recur(h, g, m, n)

%disp(['(', num2str(m), ', ', num2str(n), ')']);

M = numel(h);
N = numel(g);
MN = M * N;

if m > 1 
    if n < m
        
        f = f_recur(h, g, m - 1, n);

    else %n >= m

        f = 0;

        for k = 1:(n - m + 1)
            idx = n - k + 1;
            if k <= N
                f = f + f_recur(h, g, m - 1, idx) * g(k);
            end            
        end        
    end

else %m = 1
    
    f = 0;
    for k = 1:N
        idx = n - k + 1;
        if idx >= 1 && idx <= M
            f = f + h(idx) * g(k);
        end
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = f_recur_swap(h, g, m, n)

%disp(['(', num2str(m), ', ', num2str(n), ')']);

M = numel(h);
N = numel(g);
MN = M * N;

if m > 1 
    if n < m
        
        f = f_recur_swap(h, g, m - 1, n);

    else %n >= m

        f = 0;        

        for k = m:n
            idx = n - k + 1;
            if idx <= N
                f = f + f_recur_swap(h, g, m - 1, k) * g(idx);
            end            
        end        
    end

else %m = 1
    
    f = 0;
    for k = 1:M
        idx = n - k + 1;
        if idx >= 1 && idx <= N
            f = f + h(k) * g(idx);
        end
    end

end
