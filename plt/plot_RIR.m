function [h_fig, h_edc] = plot_RIR(h, options)
%Plot RIR, spectrogram, optional echo decay curve

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%h:     [M x 1] RIR, M taps

%options: struct, see plot_RIR_opts.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h_fig:     Handle to figure
%h_edc:     Handle to EDC curve

arguments
    h (:,1) double = [1];

    options = plot_RIR_opts;

end

h_real = real(h);
h_imag = imag(h);

if all(isreal(h)) %RIR h is all reals
    num_cols = 1;

    name_real = options.name;
    name_imag = options.name;
    
else
    num_cols = 2;

    name_real = [options.name, ' Real'];
    name_imag = [options.name, ' Imag'];
end

%Plotting
fontsize = options.font_size;

%Time-domain
h_fig = figure;
h_fig.Position = [100, 100, options.fig_size(1) * num_cols, options.fig_size(2)];
tiledlayout(options.disp_RIR + options.disp_spec, num_cols,  'TileIndexing', 'columnmajor'); 

for nc = 1:num_cols
    if nc == 1
        h = h_real;
        name = name_real;
    else
        h = h_imag;
        name = name_imag;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %RIR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if options.disp_RIR
        nexttile;
        [dst, krt, stdw] = estimate_RIR_density(h, options.Fs, options.win_size);
        edc = flipud(cumsum(flipud(h).^2)/(sum(h.^2) + eps));
        edc_dB = 10 * log10(edc);
    
        idx_p50 = find(edc <= 0.5, 1);
    
        N_h = numel(h);
        % idx_lo = ceil(N_h * options.T_frac_lo);
        % idx_hi = ceil(N_h * options.T_frac_hi);
    
        idx_hi = find(edc_dB <= options.RT60_dB_hi, 1);
        idx_lo = find(edc_dB <= options.RT60_dB_lo, 1);   
    
        t = (0:numel(h)-1)' / options.Fs;
        t_ms = t * 1000;
    
        % slope = (options.RT60_dB_hi - options.RT60_dB_lo) / (t(idx_hi) - t(idx_lo));
        % RT60 = -60 / slope;
    
        [p, stats] = polyfit(t(idx_hi:idx_lo), edc_dB(idx_hi:idx_lo), 1);
        RT60 = -60 / p(1);
            
        yyaxis left;
        plot(t_ms, h, '-', t_ms, edc, '-.', t_ms(idx_p50), 0, 'r*', 'linewidth', 1.5);
        grid on;  axis tight;
        xlabel('Time (ms)', 'fontsize', fontsize);
        ylabel('Amplitude', 'fontsize', fontsize);
        ylim([min(h) - eps, max(h)] * 1.1)
        
        yyaxis right;
        plot(t_ms, dst, 'linewidth', 1.5);
        ylim([0, max([1, max(dst)])])
        ylabel('Echo Density Profile', 'fontsize', fontsize);
        
        %Legend
        h_lg = legend('RIR', 'EDC', 'P50', 'Echo Density', 'location', options.legend_location, 'NumColumns', 2, BackgroundAlpha=.7);
        set(h_lg, 'fontsize', fontsize - 1);
        title(h_lg, ['RT60: ', num2str(RT60, 3), ' s, Rsq: ', num2str(stats.rsquared, 3)]);
        
        title([name, ' RIR'], 'fontsize', fontsize + 1);
        set(gca, 'fontsize', fontsize - 1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Spectrogram
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if options.disp_spec
        nexttile;
    
        spectrogram(h, options.win_size, ceil(options.win_size * 0.9), options.N_FFT, options.Fs, "power", "yaxis");
        axis tight; 
        set(gca, 'YScale', options.spec_scale);
        set(gca, 'fontsize', fontsize - 1);
        title([name, ' Spectrogram'], 'fontsize', fontsize + 1, 'Interpreter', options.spec_title_interp);
       
        colormap(options.colormap)
        clim(options.clim);
        if options.spec_disp_colorbar
            cb = colorbar; 
            ylabel(cb, 'Power (dB)','FontSize', fontsize);
        else
            colorbar off
        end
        xlim(options.spec_xlim);
        ylim(options.spec_ylim / 1000);

        if ~options.spec_disp_xaxis
            ax = gca;
            ax.XAxis.Visible = 'off';
        end
        if ~options.spec_disp_yaxis
            ax = gca;
            ax.YAxis.Visible = 'off';
        end

    end
end

%EDC plot
if options.disp_EDC_fig

    h_edc = figure; 
    h_edc.Position = [300, 100, options.fig_size(1) * num_cols, options.fig_size(2)];
    tiledlayout(1, num_cols,  'TileIndexing', 'columnmajor'); 


    for nc = 1:num_cols
        if nc == 1
            h = h_real;
            name = name_real;
        else
            h = h_imag;
            name = name_imag;
        end

        %RIR
        nexttile;
        [dst, krt, stdw] = estimate_RIR_density(h, options.Fs, options.win_size);
        edc = flipud(cumsum(flipud(h).^2)/sum(h.^2));
        edc_dB = 10 * log10(edc);
    
        idx_p50 = find(edc <= 0.5, 1);
    
        N_h = numel(h);
        % idx_lo = ceil(N_h * options.T_frac_lo);
        % idx_hi = ceil(N_h * options.T_frac_hi);
    
        idx_hi = find(edc_dB <= options.RT60_dB_hi, 1);
        idx_lo = find(edc_dB <= options.RT60_dB_lo, 1);   
    
        t = (0:numel(h)-1)' / options.Fs;
        
        t_ms = t * 1000;
        plot(t_ms, 10*log10(edc), t_ms(idx_lo), 10*log10(edc(idx_lo)), 'ro', t_ms(idx_hi), 10*log10(edc(idx_hi)), 'ko', 'linewidth', 1.5); grid on; axis tight;
        grid on; axis tight;
        xlabel('Time (ms)', 'fontsize', fontsize);
        ylabel('dB', 'fontsize', fontsize);
        title([name, ' Energy Decay Curve'], 'fontsize', fontsize + 1);

    end

else
    h_edc = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimate echo density of RIR over sliding window
%Echo ~ number of samples outside sliding window's std

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%h:             [N x 1] filter
%Fs:            Sampling rate
%win_size:      Window size (odd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%dst:           [N x 1] Density profile (0 = sparse, > 1 dense) over time at t
%krt:           [N x 1] Kurtosis of window at t
%stdw:          [N x 1] standard deviation of window at t
%t:             [N x 1] time per second in seconds

function [dst, krt, stdw, t] = estimate_RIR_density(h, Fs, win_size)

if nargin < 3 || isempty(win_size)
%       win_size = 64;
        win_size = 128;
%       win_size = 256;
%       win_size = 512;
        win_size = win_size + 1;
end

if rem(win_size, 2) == 0 %Force odd
    win_size = win_size + 1;
end

h = h(:);
N = numel(h);      

t = (0:N-1) / Fs;
dst = zeros(N, 1);   %Density
krt = zeros(N, 1);   %Excess Kurtosis profile
stdw = zeros(N, 1);    %Standard deviation profile

%Sliding window centered on sample n
win_size_half = (win_size - 1) / 2;
for n = 1:N    
    n_start = max(n - win_size_half, 1);
    n_end   = min(n + win_size_half, N);
    x = h(n_start:n_end);
    
%    stdw(n) = std(x);
%    dst(n) = sum((abs(x-mean(x)) > stdw(n)) / (n_end - n_start + 1) ); %Fraction of samples outside 1 std.
 
    stdw(n) = sqrt(sum(x.^2) / (n_end - n_start + 1) );
    dst(n) = sum((abs(x) > stdw(n)) / (n_end - n_start + 1) ); %Fraction of samples outside 1 std.

    krt(n) = kurtosis(x);
end
dst = dst / erfc(1/sqrt(2));  %Denominator is expected number samples lying outside 1 std

