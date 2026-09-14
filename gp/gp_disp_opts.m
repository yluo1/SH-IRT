function options = gp_disp_opts(options)
%Get default GP display options

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input and output

arguments
    options.disp_xlim (1,2) double {mustBeNonnegative} = [20, inf];
    options.disp_ylim (1,2) double {mustBeNonnegative} = [0, inf];
    options.disp_position  (1,4) double {mustBeNonnegative} =  [100, 100, 560, 480 * (0.666)];
    options.disp_legend_loc (1,:) char = 'northeast';
    options.disp_legend_num_cols (1,1) double {mustBePositive} = 1;
    options.disp_legend_samples (1,1) logical = true; %If false, display only 1 entry in legend
    options.disp_legend_mean (1,1) logical = true;  %If false, display only 1 entry in legend
    options.disp_legend_var (1,1) logical = true;   %If false, display only 1 entry in legend
    options.disp_legend_compact (1,1) logical = false; 
    options.disp_legend_transparency (1,1) double {mustBeNonnegative} = 1.0;
    options.disp_fontsize (1,1) double  {mustBePositive}  = 18;
    options.disp_mean (1,1) logical = true;
    options.disp_var (1,1) logical = true;
    options.disp_sample_eval (1,1) logical = true;
    options.disp_var_transparency (1,1) double {mustBeNonnegative} = 0.25; %variance fill transparency
    options.disp_colororder (1,:) char = 'gem'; %colororder palette
    options.disp_sample_stride (1,1) double {mustBePositive, mustBeInteger} = 1; %Subsample data when > 1
    options.disp_marker_size (1,1) double {mustBePositive} = 12; %Marker size
end