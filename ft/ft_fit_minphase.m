function [h] = ft_fit_minphase(X_mag, w, num_taps, options)
% Fit minimum phase FIR to real cepstrum of interpolated magnitude target X_mag at angular frequency w

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X_mag:         [1 x M] Target magnitude responses at angular frequency w
%w:             [1 x M] Angular frequency (radians / sample) [0 to 2 * pi]
%num_taps:      Number of FIR taps

%options:               struct
%options.rceps_M:       Number of uniform frequency points from DC to Nyquist for interpolation
%options.rceps_interp:  String, interpolation method of log-magnitude responses {'logpchip', 'loglin', 'pchip'}
%                               'logpchip'  Log-frequency transformed piece-wise cubic spline
%                               'loglin'    Log-frequency transformed piece-wise linear interp
%                               'pchip'     linear-frequency piece-wise spline
%                               'lin'       linear-frequency linear interp

%options.rceps_extrap:  String, extrapolation method at DC, Nyquist {'nearest'}
%options.rceps_fit:     String, fitting method to minimum phase target {'trunc', 'leastsquares'}
%                               'trunc'             Truncates rceps minimum phase filter
%                               'leastsquares'      Least squares re-fit to minimum phase target

%options.rceps_fit_M:   Number of log-frequency points to fit to minimum phase target 

%options.enable_disp:   Logical, if true, display fitted response
%options.Fs:            Sampling rate for display
%options.disp_dB_lim:   [1 x 2] dB limits for display [min, max]
%options.disp_deg_lim:  [1 x 2] Degree limits for display [min, max]
%options.disp_grp_del_lim:  [1 x 2] Group delay sample limits for display [min, max]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h:             [1 x num_taps]  Fitted FIR filter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Sample fit to log-frequency targets

% Fs = 48000;
% freq = [0, 100, 1000, 10000, Fs/2];
% X_mag = db2mag([2, 3, 6, -2, -4]);

% [h_16]    = ft_fit_minphase(X_mag, freq / Fs * (2 * pi), 16, 'enable_disp', true);
% [h_16_ls] = ft_fit_minphase(X_mag, freq / Fs * (2 * pi), 16, 'rceps_fit', 'leastsquares', 'enable_disp', true);

% [h_128]       = ft_fit_minphase(X_mag, freq / Fs * (2 * pi), 128, 'enable_disp', true);
% [h_128_ls]    = ft_fit_minphase(X_mag, freq / Fs * (2 * pi), 128, 'rceps_fit', 'leastsquares', 'enable_disp', true);


arguments
    X_mag (1,:) double {mustBePositive} = 1;
    w (1,:) double {mustBeNonnegative} = 0;
    num_taps (1,1) double {mustBeInteger, mustBePositive} = 1;

    options.rceps_M (1,1) double {mustBeInteger, mustBePositive} = 4096;
    options.rceps_interp (1,:) char {mustBeMember(options.rceps_interp, {'logpchip', 'loglin', 'pchip', 'lin'})} = 'logpchip';
    options.rceps_extrap (1,:) char {mustBeMember(options.rceps_extrap, {'nearest'})} = 'nearest';
    options.rceps_fit    (1,:) char {mustBeMember(options.rceps_fit, {'trunc', 'leastsquares'})} = 'trunc';
    options.rceps_fit_M (1,1) double {mustBeInteger, mustBePositive} = 4096;

    options.enable_disp (1,1) logical = false;
    options.Fs (1,1) double {mustBePositive} = 48000;
    options.disp_dB_lim (1,2) double = [-inf, inf];
    options.disp_deg_lim (1,2) double = [-180, 180];
    options.disp_grp_del_lim (1,2) double = [-10, 10];

end

assert(numel(X_mag) == numel(w), 'X_mag and w size mismatch');
assert(issorted(w), 'w must be in ascending order');
%M = numel(X_mag);

% Handle DC and Nyquist end-points for extrapolation
if strcmp(options.rceps_extrap, 'nearest')
    if w(1) ~= 0 % Add DC
        w = [0, w];
        X_mag = [X_mag(1), X_mag];
    end
    if w(end) ~= pi % Add Nyquist
        w = [w, pi];
        X_mag = [X_mag, X_mag(end)];
    end
else
    error('Unsupported options.rceps_extrap');
end

% Uniform angular frequency sampling
w_uni = linspace(0, pi, options.rceps_M);

% Handle DC as if 1 octave lower than smallest positive frequency bin
w_uni(1) = w_uni(2) / 2;
w(1) = w_uni(1);

% Interpolate
if strcmp(options.rceps_interp, 'logpchip')
    log_X_mag = interp1(log(w), log(X_mag), log(w_uni), 'pchip');
elseif strcmp(options.rceps_interp, 'loglin')
    log_X_mag = interp1(log(w), log(X_mag), log(w_uni), 'linear');  
elseif strcmp(options.rceps_interp, 'pchip')
    log_X_mag = interp1(w, log(X_mag), w_uni, 'pchip');        
elseif strcmp(options.rceps_interp, 'lin')
    log_X_mag = interp1(w, log(X_mag), w_uni, 'linear');        
else
    error('Unsupported options.rceps_interp');
end

% Compute rceps minimum phase target
xhat = real(ifft( [log_X_mag, fliplr(log_X_mag(2:end-1))] ))';
N = numel(xhat);
odd = fix(rem(N,2));
wn = [1; 2 * ones([(N+odd)/2-1,1]) ; ones([1-rem(N,2),1]); zeros([(N+odd)/2-1,1])];
H_mp_rceps = exp(fft(wn.*xhat)); % Minimum phase target
h_mp_rceps = real(ifft(H_mp_rceps))';

if strcmp(options.rceps_fit, 'trunc')

    h = h_mp_rceps(1:num_taps);
    
elseif strcmp(options.rceps_fit, 'leastsquares')

    ndx = 0:(num_taps - 1);
    w_fit = logspace(log10(w(1)), log10(pi), options.rceps_fit_M);
    F = exp( -1i * (w_fit(:) * ndx) ); %[rceps_fit_M x num_taps]
    X_tgt = freqz(h_mp_rceps, 1, w_fit);

    b = real(F'*F) \ real(F' * X_tgt(:));
    h = b';

else
    error('Unsupported options.rceps_fit');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.enable_disp && coder.target('MATLAB')

    fontsize = 16;
    linewidth = 2;
    
    hz_disp = logspace(log10(min(20, min(w(w > 0)) / (2*pi) * options.Fs)), log10(options.Fs/2), 1024);
    H_disp = freqz(h, 1, hz_disp, options.Fs);
    freq_in = w / (2*pi) * options.Fs;
    H_in = freqz(h, 1, freq_in, options.Fs);
    grp_del = grpdelay(h, 1, hz_disp, options.Fs);

    H_mp_tgt = freqz(h_mp_rceps, 1, hz_disp, options.Fs);
    grp_del_mp_tgt = grpdelay(h_mp_rceps, 1, hz_disp, options.Fs);

    h_fig = figure;
    h_fig.Position = [100,100, 900, 720];
    h_t = tiledlayout(3, 1);
    title(h_t, ['Minimum Phase Fit: ', num2str(num_taps), ' Taps'], 'fontsize', fontsize + 2);

    % Magnitude
    nexttile;
    semilogx(freq_in, mag2db(X_mag), 'ro-', 'markersize', 4, 'linewidth', linewidth); hold on;
    semilogx(hz_disp, mag2db(abs(H_disp)), 'b-', freq_in, mag2db(abs(H_in)) - mag2db(X_mag), 'g-', 'linewidth', 1.5 * linewidth); hold on;
    semilogx(hz_disp, mag2db(abs(H_mp_tgt)), 'm:', 'linewidth', 1.25 * linewidth); hold on;
    yline(0, 'k:', 'linewidth', 0.8 * linewidth);
    h_lg =  legend('Target Response', 'Fitted Response', 'Residual', 'Minimum Phase Target', ...
        'location', 'best', 'Orientation', 'horizontal');  set(h_lg, 'fontsize', fontsize - 1);
    xlabel('Frequency (Hz)', 'fontsize', fontsize); 
    ylabel('Magnitude (dB)', 'fontsize', fontsize);
    title('Magnitude Response', 'fontsize', fontsize + 1);
    grid on; axis tight;
    set(gca, 'fontsize', fontsize - 1);
    ylim(options.disp_dB_lim);
    
    % Phase
    nexttile;
    semilogx(hz_disp, rad2deg(angle(H_disp)), 'b-', 'linewidth', 1.5 * linewidth); hold on;
    semilogx(hz_disp, rad2deg(angle(H_mp_tgt)), 'm:', 'linewidth', linewidth); hold on;

    h_lg = legend('Fitted Response', 'Minimum Phase Target', 'location', 'best');  set(h_lg, 'fontsize', fontsize - 1, 'Orientation', 'horizontal');
    xlabel('Frequency (Hz)', 'fontsize', fontsize);
    ylabel('Phase (Degrees)', 'fontsize', fontsize);
    title('Phase Response', 'fontsize', fontsize + 1);
    grid on; axis tight;
    set(gca, 'fontsize', fontsize - 1);
    ylim(options.disp_deg_lim);

    % Group delay
    nexttile;
    semilogx(hz_disp, grp_del, 'linewidth', 1.5 * linewidth); hold on;
    semilogx(hz_disp, grp_del_mp_tgt, 'm:', 'linewidth', linewidth); hold on;

    h_lg = legend('Fitted Response', 'Minimum Phase Target', 'location', 'best');  set(h_lg, 'fontsize', fontsize - 1, 'Orientation', 'horizontal');
    xlabel('Frequency (Hz)', 'fontsize', fontsize);
    ylabel('Samples', 'fontsize', fontsize);
    title('Group Delay Response', 'fontsize', fontsize + 1);
    grid on; axis tight;
    set(gca, 'fontsize', fontsize - 1);
    ylim(options.disp_grp_del_lim);

end