function f = sh_int(C, mode)
%Spherical harmonic expansion's integral over unit sphere

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:         [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)

%mode:      String, integration function over spherical coordinates \Omega = (theta, phi)

%           'Sum':          Integrate SH expansion
%                               int_\Omega C.' * Y(\Omega) d\Omega
%           'SumSq':        Integrate square magnitude of expansion
%                               int_\Omega abs( C.' * Y(\Omega) )^2 d\Omega
%           'PowAvg':       Power average (normalize SumSq)
%                               int_\Omega abs( C.' * Y(\Omega) )^2 d\Omega / (4 * pi)
%           'Var':          Variance        
%                               int_\Omega abs( (C - [C(1, :), zeros((P + 1)^2 -1, M)]).' * Y(\Omega) )^2 d\Omega / (4 * pi)
%           'Entr'          Entropy (Numeric quadrature fibonnaci spiral on sphere)
%                               -int_\Omega C.' * |Y(\omega)| * log(|C.' * Y(\omega)|) d\Omega

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%f:         [1 x M] Integral

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_int', '-o', 'sh/sh_int_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare analytic versus numeric integration

% rng(648);
% M = 3;
% P = 4;
% C = sh_rand(P, M);
% [theta, phi] = sh_fib(10000);

% f_sum_int = sh_int(C, 'Sum');
% f_sum_num = mean(sh_dec(C, theta, phi), 1) * (4 * pi);
% err_sum = f_sum_int - f_sum_num

% f_sumsq_int = sh_int(C, 'SumSq');
% f_sumsq_num = mean(abs(sh_dec(C, theta, phi)).^2, 1) * (4 * pi);
% err_sumsq = f_sumsq_int - f_sumsq_num

% f_powavg_int = sh_int(C, 'PowAvg');
% f_powavg_num = mean(abs(sh_dec(C, theta, phi)).^2, 1);
% err_powavg = f_powavg_int - f_powavg_num

% f_var_int = sh_int(C, 'Var');
% f_var_num = var(sh_dec(C, theta, phi), 1);
% err_var = f_var_int - f_var_num

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compute entropy of probability density functions

% max_odr = 10;
% theta = pi/2 ;
% phi = 0;
% is_real = false;
% dB_lim = [-20 0];
% C_pdf_1 = sh_enc_rbf('SqExp', max_odr, theta, phi, 1, is_real);
% sh_plt(C_pdf_1, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Bandwidth = 1');
% f_entr_int_1 = sh_int(C_pdf, 'Entr')

% C_pdf_2 = sh_enc_rbf('SqExp', max_odr, theta, phi, 0.5, is_real);
% sh_plt(C_pdf_2, 'mercator', is_real, 'dB_lim', dB_lim, 'title_name', 'Bandwidth = 0.5');
% f_entr_int_2 = sh_int(C_pdf_2, 'Entr')

arguments
      C (:,:) double {coder.mustBeComplex} = complex(0);
      mode (1,:) char {mustBeMember(mode, {'Sum', 'SumSq', 'PowAvg', 'Var', 'Entr'})} = 'Sum';
end

if strcmp(mode, 'Sum')             %Integrate 0th order SH
    f = 2 * sqrt(pi) * C(1, :);    

elseif strcmp(mode, 'SumSq')       %Integrate over squared magnitude
    f = sum( conj(C) .* C, 1);

elseif strcmp(mode, 'PowAvg')      %Integrate over squared magnitude and normalize
    f = sum( conj(C) .* C, 1) / (4 * pi);

elseif strcmp(mode, 'Var')         %Compute variance
    f = sum( conj(C(2:end, :)) .* C(2:end, :), 1) / (4 * pi);

elseif strcmp(mode, 'Entr')        %Approximate entropy
    [theta, phi] = sh_fib(10000);
    y = abs(sh_dec(C, theta, phi));
    f = -mean(y .* log(y), 1) * (4 * pi);

else
    f = nan(1, size(C, 2));
    error('Unsupported mode');
end