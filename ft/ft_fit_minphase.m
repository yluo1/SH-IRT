function [h] = ft_fit_minphase(X_mag, w, num_taps, mode, options)
% Fit minimum phase FIR to magnitude target X_mag at angular frequency w

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X_mag:         [1 x M] Target magnitude responses at angular frequency w
%w:             [1 x M] Angular frequency (radians / sample) [0 to 2 * pi]
%num_taps:      Number of FIR taps
%mode:          String, fitting method {'rceps'}
%                   'rceps':    Fit to real cepstrum from uniform frequency interpolated X_mag

%options:               struct
%options.rceps_M:       Number of uniform frequency points from DC to Nyquist
%options.rceps_interp:  String, interpolation method {'logpchip', 'loglin', 'pchip'}
%options.rceps_extrap:  String, extrapolation method at DC, Nyquist {'nearest'}

%options.enable_disp:   Logical, if true, display fitted response
%options.Fs:            Sampling rate for display

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h:             [1 x num_taps]  Fitted FIR filter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Sample fit to log-frequency targets

% Fs = 48000;
% freq = [0, 100, 1000, 10000, Fs/2];
% X_mag = db2mag([2 3, 6, -2, -4]);
% mode = 'rceps';

% [h_16]    = ft_fit_minphase(X_mag, freq / Fs * (2 * pi), 16, mode, 'enable_disp', true);
% [h_128]   = ft_fit_minphase(X_mag, freq / Fs * (2 * pi), 128, mode, 'enable_disp', true);

arguments
    X_mag (1,:) double {mustBePositive} = 1;
    w (1,:) double {mustBeNonnegative} = 0;
    num_taps (1,1) double {mustBeInteger, mustBePositive} = 1;
    mode (1,:) char {mustBeMember(mode, {'rceps'})} = 'rceps';

    options.rceps_M (1,1) double {mustBeInteger, mustBePositive} = 4096;
    options.rceps_interp (1,:) char {mustBeMember(options.rceps_interp, {'logpchip', 'loglin', 'pchip'})} = 'logpchip';
    options.rceps_extrap (1,:) char {mustBeMember(options.rceps_extrap, {'nearest'})} = 'nearest';

    options.enable_disp (1,1) logical = false;
    options.Fs (1,1) double {mustBePositive} = 48000;
end

assert(numel(X_mag) == numel(w), 'X_mag and w size mismatch');
assert(issorted(w), 'w must be in ascending order');

if strcmp(mode, 'rceps') % Real cepstrum method
    
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

    %Interpolate
    if strcmp(options.rceps_interp, 'logpchip')
        log_X_mag = interp1(log(w), log(X_mag), log(w_uni), 'pchip');
    elseif strcmp(options.rceps_interp, 'loglin')
        log_X_mag = interp1(log(w), log(X_mag), log(w_uni), 'linear');  
    elseif strcmp(options.rceps_interp, 'pchip')
        log_X_mag = interp1(w, log(X_mag), w_uni, 'pchip');        
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
    h = h_mp_rceps(1:num_taps);

else
    error('Unsupported mode');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.enable_disp && coder.target('MATLAB')

    fontsize = 16;
    
    hz_disp = logspace(log10(min(20, w(1) / (2*pi) * options.Fs)), log10(options.Fs/2), 1024);
    H_disp = freqz(h, 1, hz_disp, options.Fs);
    grp_del = grpdelay(h, 1, hz_disp, options.Fs);

    h_fig = figure;
    h_fig.Position = [100,100, 800, 600];
    h_t = tiledlayout(3, 1);
    title(h_t, ['Minimum Phase Fit: ', num2str(num_taps), ' Taps'], 'fontsize', fontsize + 2);

    % Magnitude
    nexttile;
    semilogx(hz_disp, mag2db(abs(H_disp)), w / (2*pi) * options.Fs, mag2db(X_mag), '*', 'markersize', 12, 'linewidth', 1.5);
    h_lg =  legend('Fitted Response', 'Target Response', 'location', 'best');  set(h_lg, 'fontsize', fontsize - 1);
    xlabel('Frequency (Hz)', 'fontsize', fontsize); 
    ylabel('Magnitude (dB)', 'fontsize', fontsize);
    title('Magnitude Response', 'fontsize', fontsize + 1);
    grid on; axis tight;
    set(gca, 'fontsize', fontsize - 1);
    
    % Phase
    nexttile;
    semilogx(hz_disp, deg2rad(angle(H_disp)), 'linewidth', 1.5);
    h_lg = legend('Fitted Response', 'location', 'best');  set(h_lg, 'fontsize', fontsize - 1);
    xlabel('Frequency (Hz)', 'fontsize', fontsize);
    ylabel('Phase (Degrees)', 'fontsize', fontsize);
    title('Phase Response', 'fontsize', fontsize + 1);
    grid on; axis tight;
    set(gca, 'fontsize', fontsize - 1);

    % Group delay
    nexttile;
    semilogx(hz_disp, grp_del, 'linewidth', 1.5);
    h_lg = legend('Fitted Response', 'location', 'best');  set(h_lg, 'fontsize', fontsize - 1);
    xlabel('Frequency (Hz)', 'fontsize', fontsize);
    ylabel('Samples', 'fontsize', fontsize);
    title('Group Delay Response', 'fontsize', fontsize + 1);
    grid on; axis tight;
    set(gca, 'fontsize', fontsize - 1);

end