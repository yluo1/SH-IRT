function options = plot_RIR_opts(options)
%Get plot_RIR options

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input / output

arguments

    options.Fs (1,1) double {mustBePositive, mustBeInteger} = 48000;
    options.win_size (1,1) double {mustBePositive, mustBeInteger}= 4096;
    options.N_FFT (1,1) double {mustBePositive, mustBeInteger}= 4096;

    options.fig_size (1,2) double {mustBePositive} = [900 600] * (3/4);
    options.font_size (1,1) double {mustBePositive} = 16;

    options.name (1,:) char = '';

    options.spec_scale (1,:) char {mustBeMember(options.spec_scale, {'linear', 'log'})}  = 'log';
    options.spec_disp_colorbar (1,1) logical = true;
    options.spec_disp_xaxis (1,1) logical = true;
    options.spec_disp_yaxis (1,1) logical = true;
    options.spec_title_interp (1,:) char {mustBeMember(options.spec_title_interp, {'none', 'latex'})} = 'none';
    options.colormap = parula;
    options.clim (1,2) double = [-80, -20];
    options.spec_ylim (1,2) double = [-inf inf];
    options.spec_xlim (1,2) double = [-inf inf];
    options.spec_xunits (1,:) char {mustBeMember(options.spec_xunits, {'ms', 'sec'})} = 'ms'

    options.legend_location (1,:) char = 'east';

    options.RT60_dB_hi (1,1) double = -10;
    options.RT60_dB_lo (1,1) double = -30;    

    options.disp_RIR (1,1) logical = true;
    options.disp_spec  (1,1) logical = true;
    options.disp_EDC_fig (1,1) logical = false;

end