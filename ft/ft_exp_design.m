function [g, g_minphase, h_f] = ft_exp_design(num_taps, mode, options)
%Design exponentiating FIR filter g

%Author: Yuancheng Luo, 2026

%Paper Reference:
%Yuancheng Luo, "Fast Time-Varying Exponentiated Convolution Methods for Generative Direction Dependent Reverberation",
%Proceedings of the 161th Audio Engineering Society Convention.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%mode:          String, specify filter design targets {'RT60', 'air_absorption_bass', 'air_absorption_iso'}
%                   'RT60':                 Room attenuation seconds by 60 dB in options.RT60_sec
%                   'air_absorption_iso'    Air-absorption attenuation (dB / m) ISO 9613-1:1993
%                   'air_absorption_bass'   Air-absorption attenuation (dB / m) Bass, H. E. Reference

%num_taps:      Number of taps in FIR filter design

%options:       struct

%options.fit_method:    String, filter fit method {'constr_min_phase', 'constr_min_phase_minimax', 'two_tap'}
%                           'constr_min_phase_ls':       Magnitude upper-bound constrained minimum phase least squares
%                           'constr_min_phase_minimax':  Magnitude upper-bound constrained minimum phase minimax
%                           'two_tap':                   DC and Nyquist constrained two-tap FIR

%options.Fs:    Sampling rate

%options.air_params:                    struct
%options.air_params.T:                  Temperature (Celsius)            
%options.air_params.ps:                 Relative Humidity (Percent)
%options.air_params.num_bins_DC_NQ:     Number of uniformed spaced frequency bins between DC and Nyquist

%options.RT60_sec:      [1 x num_bins_DC_NQ] time (seconds) for 60 dB attenuation at uniform spaced frequency bin centers between DC and Nyquist

%options.tol0:          Margin of exclusion below unity for filter's magnitude response

%options.enable_disp:   Logical, if true, plot filter responses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%g:             [1 x num_taps] FIR
%g_minphase:    [1 x *] FIR minphase target filter
%h_f:           Handle to figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Atmosphere attenuation

%g_bass = ft_exp_design(6, 'air_absorption_bass', 'enable_disp', true, 'Fs', 96000);
%g_iso  = ft_exp_design(6, 'air_absorption_iso', 'enable_disp', true, 'Fs', 96000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Custom RT60
%g_RT60 = ft_exp_design(6, 'RT60',  'enable_disp', true, 'Fs', 48000, 'RT60_sec', [5 0.35 0.9 0.2 0.3]);

%With overshoot
%g_RT60 = ft_exp_design(6, 'RT60',  'enable_disp', true, 'Fs', 48000, 'RT60_sec', [5 0.35 8 4 0.3]);

arguments
    num_taps double {mustBeInteger, mustBePositive(num_taps)} = 6;       %Number of output filter taps

    mode (1,:) char {mustBeMember(mode, {'air_absorption_bass', 'air_absorption_iso', 'RT60'})} =  'RT60';

    options.fit_method (1,:) char {mustBeMember(options.fit_method, {'constr_min_phase_ls', 'constr_min_phase_minimax', 'two_tap'}) } = 'constr_min_phase_ls';
    options.Fs {mustBePositive(options.Fs)} = 48000;     
    options.air_params = struct('T', 25, 'hr', 50, 'ps', 1, 'num_bins_DC_NQ', 129);
    options.RT60_sec (1,:) double {mustBePositive} =  [1 0.5];
    options.tol0 (1,1) double {mustBeNonnegative} = 0; 
    options.enable_disp (1,1) logical = false;
end

%Specify filter target:
%H_tgt_dB_DC_NQ
%H_tgt_abs_DC_NQ

meters_per_sample = 343 / options.Fs;
if contains(mode, 'air_absorption') %air_absorption_bass, air_absorption_iso
    
    M_half = options.air_params.num_bins_DC_NQ; %Number of bins from DC to Nyquist    
    freq = linspace(0, options.Fs / 2, M_half); %DC to Nyquist, uniform spaced frequency

    %Compute dB attenuation / meter
    dB_bass_m_Hz = zeros(1, M_half);
    dB_iso_m_Hz = zeros(1, M_half);
    for i = 1:M_half
        [dB_bass_m_Hz(i), dB_iso_m_Hz(i)] = air_absorption(freq(i), options.air_params.T,  options.air_params.hr,  options.air_params.ps);
    end
    if strcmp(mode, 'air_absorption_bass')
        dB_m_Hz = -dB_bass_m_Hz;
    else
        dB_m_Hz = -dB_iso_m_Hz;
    end

    H_tgt_dB_DC_NQ  = dB_m_Hz * meters_per_sample;  %Target magnitude dB
    H_tgt_abs_DC_NQ = db2mag(H_tgt_dB_DC_NQ);       %Target magnitude modulus    

elseif strcmp(mode, 'RT60') %Compute filter to achieve target RT60

    H_tgt_dB_DC_NQ  = -60 ./ (options.Fs * options.RT60_sec); %Target magnitude dB
    H_tgt_abs_DC_NQ = db2mag(H_tgt_dB_DC_NQ);                 %Target magnitude modulus  

else
    error('Unknown mode')
end


% Fit filter
if contains(options.fit_method, 'constr_min_phase')
    X_mag = [H_tgt_abs_DC_NQ, fliplr(H_tgt_abs_DC_NQ(2:end-1))];

    if strcmp(options.fit_method, 'constr_min_phase_ls')
        mode_ft_bnd_minphase = 'least_squares';
    elseif strcmp(options.fit_method, 'constr_min_phase_minimax')
        mode_ft_bnd_minphase = 'minimax';
    else
        error('Unsupported mode');
    end
    [g, g_minphase, err, lambda] = ft_bnd_minphase(X_mag(:), min(numel(X_mag), num_taps), 'ub', 1 - options.tol0, 'mode', mode_ft_bnd_minphase);

    g = g(:)';
    g = [g, zeros(1, num_taps - numel(g))];

elseif strcmp(options.fit_method, 'two_tap')

    g = ft_two_tap_FIR(H_tgt_dB_DC_NQ(1), H_tgt_dB_DC_NQ(end));
    g = [g, zeros(1, num_taps - 2)];
    g = g(1:num_taps);

    g_minphase = g;

else
    error('Unsupported options.fit_method');
end

%Check margins and scale g if in exclusion region
G_margin_check = freqz(g, 1, linspace(0, options.Fs/2, 512), options.Fs);
G_margin_check_abs = abs(G_margin_check);
max_G_margin_check_abs = max(G_margin_check_abs);
if max_G_margin_check_abs > (1 - options.tol0)    
    g = g / max_G_margin_check_abs;
    warning('attenuating g as |G(w)| exceeded 1 - tol0');
end

%Plot target and filter responses
if options.enable_disp

    freq_disp = logspace(log10(20), log10(options.Fs/2), 512);
    
    G_disp = freqz(g, 1, freq_disp, options.Fs);
    H_disp = freqz(g_minphase, 1, freq_disp, options.Fs);
%   H_disp = db2mag(interp1(freq, H_tgt_dB_DC_NQ, freq_disp, 'pchip'));

    h_f = plot_cpx_resp([H_disp(:), G_disp(:)], freq_disp(:), 'resp_names', {'Target', 'Filter'});
else
    h_f = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha, alpha_iso, c, c_iso] = air_absorption(f, T, hr, ps)
% % air_absorption: Calculates sound absorption (attenuation) in humid air   
% % 
% % Syntax:
% % 
% % [alpha, alpha_iso, c, c_iso]=air_absorption(f, T, hr, ps);
% %
% % **********************************************************************
% % 
% % Description
% % 
% % Sound absorption (attenuation) in humid air depends on frequency, 
% % temperature, relative humidity, and atmospheric pressure.  
% % 
% % Calculates sound absorption (attenuation) in humid air using the ISO
% % standard and the Bass forumla.  The sound absorption is output in
% % dB/meter.
% % 
% % See appropriate input and output variables sections below for
% % more details. 
% % 
% % 
% % **********************************************************************
% %
% % Input Variables
% % 
% % f=100;      % (Hz) frequency of pure tone
% %             % default is f=100;
% %
% % T=25;       % (degrees Celsius) temperature
% %             % default is T=25;
% % 
% % hr=50;      % (Percent) Relative Humidity 
% %             % default is hr=50;
% % 
% % ps=1;       % atmospheric pressure ratio
% %             % pa/pr ratio of ambient atmospoheric pressure to the
% %             % standard atmosphere.  
% %             % default is ps=1;
% % 
% % **********************************************************************
% %
% % Output Variables
% % 
% % alpha       sound absorption is (dB/meter) according to Bass
% % 
% % alpha_iso   sound absorption is (dB/meter) according to ISO standard
% % 
% % c           speed of sound in humid air (meters/second) according to
% %             Bass
% % 
% % c_iso       speed of sound in humid air (meters/second) according to 
% %             ISO standard
% % 
% % **********************************************************************
% 
% Example='1';
% 
% 
% f=100;   % frequency  in Hz
%
% T=20;    % 20 Degrees Celsius
%
% hr=80;   % Relative humidity in percentage hr=80 means 80 percent humidity
%          
% ps=1;    % Is the barometric pressure ratio. Usually, ps=1;
% 
% % Run the program
% [alpha, alpha_iso, c, c_iso]=air_absorption(f, T, hr, ps);
% 
% 
% % **********************************************************************
% % 
% % References
% % 
% % ANSI/ASA S1.26-1995 (R2009)
% % 
% % ISO 9613-1:1993 Acoustics -- Attenuation of sound during propagation 
% % outdoors -- Part 1: Calculation of the absorption of sound by 
% % the atmosphere
% % 
% % 
% % Need both Articles below inclduing the (Further Developments) article 
% % from Bass to make the calculations correctly.  
% % 
% % Attenborough, K., S. Taherzadeh, H. E. Bass, X. Di, R. Raspet, 
% % G. R. Becker, A. Güdesen et al. "Benchmark cases for outdoor sound 
% % propagation models." The Journal of the Acoustical Society of America 
% % 97 (1995): 173.
% % 
% % Bass, H. E., L. C. Sutherland, A. J. Zuckerwar, D. T. Blackstock, 
% % and D. M. Hester. "Atmospheric absorption of sound: Further 
% % developments." The Journal of the Acoustical Society of America 
% % 97, no. 1 (1995): 680-683.
% % 
% % Pierce A.D., Acoustics an introduciton to its physical principles and
% % aplications, Acoustical Society of America, Equation 1-9.5 pp.30
% %
% %
% %
% % **********************************************************************
% %
% % This program was  Written by Edward L. Zechmann   
% % 
% % date        September   2006  
% % 
% % modified    March       2010    Updated Comments chnged input units of
% %                                 T from Farhenheit to Celsius.
% % 
% % modified    October 24  2013    Updated formula for speed of sound by 
% %                                 dividing, h the humidity in percent 
% %                                 molecular concentration by 100 to 
% %                                 convert it to a fraction.  
% % 
% % 
% % **********************************************************************
% % 
% %
% % Please Feel Free to Modify This Program
% %
% % See Also: Atmosphere (its on the Matlab file exchange)
% %
if nargin < 1 || isempty(f) || ~isnumeric(f)
    f=100;
end
if nargin < 2 || isempty(T) || ~isnumeric(T)
    T=25;
end
if nargin < 3 || isempty(hr) || ~isnumeric(hr)
    hr=50;
end
if nargin < 4 || isempty(ps) || ~isnumeric(ps)
    ps=1;
end
% convert T from Celsius to Kelvin;
%T=273.15+5/9*(T-32); to convert from Fahrenheit
T=273.15+T; % To convert from Celsius
% listing of constants
T01=273.16; %triple point in degrees Kelvin
T0=293.15;
% atmospheric pressure ratio is the ambient pressure/standard pressure
ps0=1;  % ps0= standard pressure/standard pressure which is unity
% Bass formula for saturation pressure ratio
psat_ps0=10^( 10.79586*(1-T01/T) -5.02808*log10(T/T01) +1.50474*10^(-4)*(1-10^(-8.29692*(T/T01-1))) -4.2873*10^(-4)*(1-10^(-4.76955*(T01/T-1)))-2.2195983);
% Iso formula for saturation pressure ratio
psat_ps0_iso=10^( -6.8346*(T01/T)^1.261+4.6151);
ps_ps0=ps/ps0;
h = hr*psat_ps0/ps_ps0; % h is the humidity in percent molar concentration
h_iso = hr*psat_ps0_iso/ps_ps0; % h is the humidity in percent molar concentration
c0 = 331; % c0 is the reference sound speed
c = (1+0.16*h/100)*c0*sqrt(T/T01);
c_iso = (1+0.16*h_iso/100)*c0*sqrt(T/T01);
% % 
% % **********************************************************************
%  
% 
% Bass formula
F=f/ps;
Fr0=1/ps0*(24+4.04*10^4*h*(0.02+h)/(0.391+h));
FrN=1/ps0*(T0/T)^(1/2)*(9+280*h*exp(-4.17*((T0/T)^(1/3)-1)));
% Calculate the air absorption in dB/meter using the Bass formula
alpha=20*log10(exp(1))*ps*F^2*( 1.84*10^(-11)*(T/T0)^(0.5)*ps0+((T/T0)^(-5/2))*( 0.01275*exp(-2239.1/T)/(Fr0+F^2/Fr0)+0.1068*exp(-3352/T)/(FrN+F^2/FrN) ));
% % 
% % **********************************************************************
%  
% 
%ISO formula
taur=T/T0;
pr=ps/ps0;
fr0=pr*(24+40400*h_iso*(0.02+h_iso)/(0.391+h_iso));
frN=pr*(taur)^(-1/2)*(9+280*h_iso*exp(-4.17*((taur)^(-1/3)-1)));
b1=0.1068*exp(-3352/T)/(frN+f^2/frN);
b2=0.01275*exp(-2239.1/T)/(fr0+f^2/fr0);
% Calculate the air absorption in dB/meter for the ISO standard
alpha_iso=8.686*f^2*taur^(1/2)*(1.84*10^(-11)/pr+taur^(-3)*(b1+b2));

