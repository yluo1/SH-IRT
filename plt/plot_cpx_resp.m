function h_f = plot_cpx_resp(H, hz, options)
%Plot complex frequency responses's magnitude, phase, with optional 
%impulse and group delay responses

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%H:         [M x N]          M bins x N responses
%hz:        [M x 1]          Frequency vector (Hz), non-negative

%options:       struct

%options.title_name:         String, plot title name
%options.resp_names:         [1 x N]    cell of strings or chars for names of each response, {} for default numbering

%options.show_mag:           Logical, if true, display magnitude response
%options.show_ph:            Logical, if true, display phase response
%options.show_IR:            Logical, if true, display impulse response if non-empty
%options.show_GD:            Logical, if true, display sample group delay if non-empty

%options.dB_lim:             [1 x 2]    Min and Max dB

%options.IR:                 [T x N]    Impulse response, [] to ignore 
%options.IR_x_lim:           [1 x 2]    Min and Max sample range for IR

%options.GD:                 [M x N]    Sample group delay, [] to ignore

%options.linewidth:          Scalar, line width, must be positive
%options.fontsize:           Scalar, fontsize, must be positive

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h_f:       Figure handle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Plot butterworth filter responses

% M = 512;
% N = 2;
% Fs = 48000;
% [b1,a1] = butter(8, pi/8);
% [b2,a2] = butter(8, pi/16);
% hz = logspace(log10(20), log10(Fs / 2), M)';
% H = [freqz(b1, a1, hz, Fs), freqz(b2, a2, hz, Fs)];
% plot_cpx_resp(H, hz, 'dB_lim', [-60, 1]);


arguments
    H (:,:) double = 1;
    hz (:,1) double {mustBeNonnegative} = 0;

    options.title_name (1,:) char = 'Frequency Response';
    options.resp_names (1,:) cell = {};

    options.show_mag (1,1) logical   = true;
    options.show_phase (1,1) logical = true;
    options.show_IR  (1,1) logical   = true;
    options.show_GD   (1,1) logical  = true;

    options.dB_lim  (1,2) double = [-inf inf];

    options.IR (:,:) double = [];
    options.IR_x_lim (1,2) double = [-inf inf];

    options.GD (:,:) double  = [];

    options.legend_num_rows (1,1) double {mustBePositive} = 1;

    options.linewidth (1,1) double {mustBePositive}  = 1.5;

    options.fontsize (1,1) double {mustBePositive}  = 14;

end

hz = hz(:);

[M, N] = size(H);
assert(M == numel(hz), 'H, hz size mismatch');

H_mag_dB = mag2db(abs(H));
H_deg    = rad2deg(angle(H));

%Compute number of figure tiles
num_tiles = options.show_mag + options.show_phase + ...
            (options.show_IR && ~isempty(options.IR)) + ...
            (options.show_GD && ~isempty(options.GD));

h_f = figure;
h_f.Position = [100, 100, 1024, 768/2 * num_tiles];

fontsize = options.fontsize;
     
num_cols = ceil(N / options.legend_num_rows);

h_t = tiledlayout(num_tiles, 1);
h_t.Title.String = options.title_name;
h_t.Title.FontSize = fontsize + 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.show_mag
    nexttile;
    
    semilogx(hz, H_mag_dB, 'linewidth', options.linewidth);
    grid on; axis tight;
    xlabel('Frequency (Hz)', 'fontsize', fontsize);
    ylabel('Magnitude (dB)', 'fontsize', fontsize);
    
    set(gca, 'fontsize', fontsize - 1);

    h_lg = legend(options.resp_names, 'location', 'southoutside', 'orientation', 'horizontal', 'NumColumns', num_cols);
    set(h_lg, 'fontsize', fontsize - 1);

    ylim(options.dB_lim)

end

if options.show_phase
    nexttile;

    semilogx(hz, H_deg, 'linewidth', options.linewidth);
    grid on; axis tight;
    xlabel('Frequency (Hz)', 'fontsize', fontsize);
    ylabel('Phase (Degree)', 'fontsize', fontsize);
    
    set(gca, 'fontsize', fontsize - 1);

    h_lg = legend(options.resp_names, 'location', 'southoutside', 'orientation', 'horizontal',  'NumColumns', num_cols);
    set(h_lg, 'fontsize', fontsize - 1);

end

if options.show_GD && ~isempty(options.GD)

    nexttile;

    semilogx(hz, options.GD, 'linewidth', options.linewidth);
    grid on; axis tight;
    xlabel('Frequency (Hz)', 'fontsize', fontsize);
    ylabel('Group Delay (Samples)', 'fontsize', fontsize);
    
    set(gca, 'fontsize', fontsize - 1);

    h_lg = legend(options.resp_names, 'location', 'southoutside', 'orientation', 'horizontal', 'NumColumns', num_cols);
    set(h_lg, 'fontsize', fontsize - 1);

end

if options.show_IR && ~isempty(options.IR)
    nexttile;

    T = size(options.IR, 1);
    plot(1:T, options.IR, 'linewidth', options.linewidth); 
    grid on; axis tight;
    xlabel('Time (Samples)', 'fontsize', fontsize);
    ylabel('Amplitude', 'fontsize', fontsize);

    set(gca, 'fontsize', fontsize - 1);

    h_lg = legend(options.resp_names, 'location', 'southoutside', 'orientation', 'horizontal', 'NumColumns', num_cols);
    set(h_lg, 'fontsize', fontsize - 1);

    xlim(options.IR_x_lim)

end

