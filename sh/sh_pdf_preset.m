function C_pdf = sh_pdf_preset(dir_name, func_name, max_odr, is_real, options)
%Generate spherical harmonic probability density function presets

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%dir_name:      String, direction name
%                   {'FwdCenter', 'FwdLeft', 'FwdRight', 'FwdTop', 'FwdBot', 'Behind', ...
%                       'North', 'South', 'East', 'West', ...
%                       'SouthWest', 'NorthWest', 'SouthEast', 'NorthEast'} 

%func_name:     String, function mode
%                  {'Dirac', 'SqExp', 'SqExpLeftRight', 'SqExpTopBot', 'SqExpNorthSouth', 'SqExpSlash', 'SqExpBackSlash', 'DiracRandom', 'Uniform'}

%max_odr:       Max SH order 
%is_real:       Logical, if true, C is real, otherwise, C is complex
%enable_disp:   Logical, if true, plot preset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C              [(max_odr + 1)^2 x num_func] SH coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate sample PDF presets

%max_odr = 6;
%enable_disp = true;
%is_real = false;

%sh_pdf_preset('FwdCenter', 'Dirac', max_odr, is_real, 'enable_disp', enable_disp);
%sh_pdf_preset('FwdLeft', 'Dirac', max_odr, is_real, 'enable_disp', enable_disp);
%sh_pdf_preset('FwdRight', 'Dirac', max_odr, is_real, 'enable_disp', enable_disp);
%sh_pdf_preset('FwdTop', 'Dirac', max_odr, is_real, 'enable_disp', enable_disp);
%sh_pdf_preset('FwdBot', 'Dirac', max_odr, is_real, 'enable_disp', enable_disp);

%sh_pdf_preset('FwdBot', 'SqExp', max_odr, is_real, 'enable_disp', enable_disp);
%sh_pdf_preset('FwdBot', 'SqExpLeftRight', max_odr, is_real, 'enable_disp', enable_disp);

%sh_pdf_preset('FwdBot', 'SqExpDiagLeft', max_odr, is_real, 'enable_disp', enable_disp);
%sh_pdf_preset('FwdBot', 'SqExpDiagRight', max_odr, is_real, 'enable_disp', enable_disp);


%sh_pdf_preset('North', 'SqExp', max_odr, is_real, 'enable_disp', enable_disp);
%sh_pdf_preset('South', 'SqExp', max_odr, is_real, 'enable_disp', enable_disp);
%sh_pdf_preset('East', 'SqExp', max_odr, is_real, 'enable_disp', enable_disp);
%sh_pdf_preset('West', 'SqExp', max_odr, is_real, 'enable_disp', enable_disp);

% sh_pdf_preset('SouthWest', 'SqExp', max_odr, is_real, 'enable_disp', enable_disp);
% sh_pdf_preset('NorthWest', 'SqExp', max_odr, is_real, 'enable_disp', enable_disp);
% sh_pdf_preset('SouthEast', 'SqExp', max_odr, is_real, 'enable_disp', enable_disp);
% sh_pdf_preset('NorthEast', 'SqExp', max_odr, is_real, 'enable_disp', enable_disp);

%rng(123);
%sh_pdf_preset('FwdCenter', 'DiracRandom', max_odr, is_real, 'enable_disp', enable_disp);

arguments
    dir_name (1,:) char {mustBeMember(dir_name, ...
         {'FwdCenter', 'FwdLeft', 'FwdRight', 'FwdTop', 'FwdBot', 'Behind', ...
         'North', 'South', 'East', 'West', ...
         'SouthWest', 'NorthWest', 'SouthEast', 'NorthEast'} )} = 'FwdCenter';
   
    func_name (1,:) char {mustBeMember(func_name, ...
          {'Dirac', 'SqExp', 'SqExpLeftRight', 'SqExpTopBot', 'SqExpNorthSouth', 'SqExpDiagLeft', 'SqExpDiagRight', 'DiracRandom', 'Uniform'})} = 'Dirac';

    max_odr (1,1) double {mustBeNonnegative} = 6;
    is_real (1,1) logical = false;
    
    options.ell (1,1) double {mustBePositive} = 3/4;
    options.enable_disp (1,1) logical = false;
end

max_odr_half = floor(max_odr / 2);

%Direction
if strcmp(dir_name, 'FwdCenter')
    theta = pi/2;
    phi = 0;

elseif strcmp(dir_name, 'FwdLeft')
    theta = pi/2;
    phi = -pi/3;

elseif strcmp(dir_name, 'FwdRight')
    theta = pi/2;
    phi = pi/3;

elseif strcmp(dir_name, 'FwdTop')
    theta = pi/2 - pi/3;
    phi = 0;

elseif strcmp(dir_name, 'FwdBot')
    theta = pi/2 + pi/3;
    phi = 0;

elseif  strcmp(dir_name, 'Behind')
    theta = pi/2;
    phi = pi;

elseif  strcmp(dir_name, 'North')
    theta = 0;
    phi = 0;

elseif  strcmp(dir_name, 'South')
    theta = pi;
    phi = 0;

elseif  strcmp(dir_name, 'East')
    theta = pi/2;
    phi = pi/2;

elseif  strcmp(dir_name, 'West')
    theta = pi/2;
    phi = -pi/2;

elseif  strcmp(dir_name, 'SouthWest')
    theta = pi/2 + pi/3;
    phi = -pi/3;

elseif  strcmp(dir_name, 'NorthWest')
    theta = pi/2 - pi/3;
    phi = -pi/3;

elseif  strcmp(dir_name, 'SouthEast')
    theta = pi/2 + pi/3;
    phi = pi/3;

elseif  strcmp(dir_name, 'NorthEast')
    theta = pi/2 - pi/3;
    phi = pi/3;


else
    error('Unknown dir_name');
end

%Function
ell = options.ell;
if strcmp(func_name, 'Dirac')
    C_pdf = sh_msq(sh_enc_proj(max_odr_half, theta, phi, is_real), is_real);    

elseif strcmp(func_name, 'SqExp')
    C_pdf = sh_msq(sh_enc_rbf('SqExp', max_odr_half, theta, phi, ell, is_real), is_real);    

elseif strcmp(func_name, 'SqExpLeftRight')
    C_pdf = sh_msq(sh_enc_rbf('SqExp', max_odr_half, pi/2,  pi/2, ell, is_real), is_real) ...
         +  sh_msq(sh_enc_rbf('SqExp', max_odr_half, pi/2, -pi/2, ell, is_real), is_real);

elseif strcmp(func_name, 'SqExpTopBot')
    C_pdf = sh_msq(sh_enc_rbf('SqExp', max_odr_half, pi/2 - pi/3,  0, ell, is_real), is_real) ...
         +  sh_msq(sh_enc_rbf('SqExp', max_odr_half, pi/2 + pi/3,  0, ell, is_real), is_real);

elseif strcmp(func_name, 'SqExpNorthSouth')
    C_pdf = sh_msq(sh_enc_rbf('SqExp', max_odr_half, 0,  0, ell, is_real), is_real) ...
         +  sh_msq(sh_enc_rbf('SqExp', max_odr_half, pi,  0, ell, is_real), is_real);

elseif strcmp(func_name, 'SqExpDiagLeft')
    C_pdf = sh_msq(sh_enc_rbf('SqExp', max_odr_half, pi/2 - pi/3,  pi/3, ell, is_real), is_real) ...
         +  sh_msq(sh_enc_rbf('SqExp', max_odr_half, pi/2 + pi/3,  -pi/3, ell, is_real), is_real);

elseif strcmp(func_name, 'SqExpDiagRight')
    C_pdf = sh_msq(sh_enc_rbf('SqExp', max_odr_half, pi/2 - pi/3,  -pi/3, ell, is_real), is_real) ...
         +  sh_msq(sh_enc_rbf('SqExp', max_odr_half, pi/2 + pi/3,  pi/3, ell, is_real), is_real);

elseif strcmp(func_name, 'DiracRandom')
    num_func = 20;
    is_mag_sq = true;
    C_pdf = sum(sh_rand_proj(max_odr, num_func, is_real, is_mag_sq), 2);

elseif strcmp(func_name,'Uniform')
    C_pdf = sh_enc_uni(max_odr);

else
    error('Unknown func_name');
end

%Normalize
C_pdf = sh_nrm(C_pdf, 'Sum');
C_pdf = [C_pdf; zeros((max_odr + 1)^2  - size(C_pdf, 1), 1)];

%Plotting
if options.enable_disp && coder.target('MATLAB')
    sh_plt(C_pdf, 'mercator', is_real, 'title_name', [dir_name, ' ', func_name], 'dB_lim', [-48, 6], 'disp_phase', false);
end