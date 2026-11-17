function [wt, h_fig] = ft_freq_wt_chs(freq, range, dBOct, options)
% Generate multichannel operating frequency weights

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%freq:      [1 x M] Frequency (Hz)
%range:     [N x 2] N channels of [fc_min, fc_max]
%dBOct:     [N x 2] N channels of [dBOct_min, dBOct_max]

% options:                  struct

% options.norm_sum:         Logical, if true, normalize weights to sum to unity per frequency
% options.soft_knee:        Logical, if true, dB/Oct weights have soft knees at -6 dB at fc

% options.enable_disp:      Logical, if true, display plot
% options.disp_dB_lim:      [1 x 2]     display dB range [min, max]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%wt:        [N x M] Weights
%h_fig:     Handle to figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate sample weights

%M = 256;
%freq = logspace(log10(20), log10(24000), M);
%disp_dB_lim = [-48, 0];
% wt1 = ft_freq_wt_chs(freq, [100, 10000; 0 1000; 4000 8000], [24 24; 24 24; 24 0], ... 
%                   'enable_disp', true, 'disp_dB_lim', disp_dB_lim, 'soft_knee', true, 'norm_sum', false);
% wt2 = ft_freq_wt_chs(freq, [100, 10000; 0 1000; 4000 8000], [24 24; 24 24; 24 0], ... 
%                  'enable_disp', true, 'disp_dB_lim', disp_dB_lim, 'soft_knee', true, 'norm_sum', true);
% wt3 = ft_freq_wt_chs(freq, [100, 10000; 0 1000; 4000 8000], [24 24; 24 24; 24 0], ... 
%                  'enable_disp', true, 'disp_dB_lim', disp_dB_lim, 'soft_knee', false, 'norm_sum', true);

arguments
    freq (1,:) double = [0];
    range (:,2) double {mustBeNonnegative} = [];
    dBOct (:,2) double {mustBeNonnegative} = [];
    
    options.norm_sum (1,1) logical = false;
    options.soft_knee (1,1) logical = true;

    options.enable_disp (1,1) logical = false;
    options.disp_dB_lim (1,2) double = [-inf inf];
    
end

M = numel(freq);
N = size(range, 1);
assert(isequal(size(range), size(dBOct) ), 'dBOct size mismatch range')

% Set filter types
if options.soft_knee
    mode_HP = 'octaveSoftHP';
    mode_LP = 'octaveSoftLP';
else
    mode_HP = 'octaveHP';
    mode_LP = 'octaveLP';
end

% Generate weights per channel
wt = zeros([N, M]);
for n = 1:N
    
    assert(range(n, 1) <= range(n,2), 'range(n,:) must be in ascending order');

    if range(n,1) == 0 || dBOct(n,1) == 0 % Lower half
        wt_lo_n = ft_freq_wt(freq, 'unity');
    else
        wt_lo_n = ft_freq_wt(freq, mode_HP, 'oct_freq', range(n,1), 'oct_roll', dBOct(n,1));
    end

    if range(n,2) == inf || dBOct(n,2) == 0 % Upper half
        wt_hi_n = ft_freq_wt(freq, 'unity');
    else
        wt_hi_n = ft_freq_wt(freq, mode_LP, 'oct_freq', range(n,2), 'oct_roll', dBOct(n,2));
    end

    wt(n,:) = (wt_lo_n .* wt_hi_n)';
end

% Normaize summation to unity
if options.norm_sum
    wt = bsxfun(@rdivide, wt, sum(wt, 1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_fig = [];
if options.enable_disp && coder.target("MATLAB")
    fontsize = 16;
    h_fig = figure;
    for n = 1:N
        semilogx(freq, mag2db(wt(n,:)), 'linewidth', 2 ); hold on;
    end
    grid on; axis tight;
    xlabel('Frequency (Hz)', 'fontsize', fontsize);
    ylabel('Magnitude (dB)', 'fontsize', fontsize);
    title('Channel Frequency Weights', 'fontsize', fontsize + 1); 
    set(gca, 'fontsize', fontsize - 1);
    ylim(options.disp_dB_lim);
end

