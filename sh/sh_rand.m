function C = sh_rand(max_odr, num_func, is_real, is_real_dec)
%Generate normally distributed random function over sphere

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%max_odr:       Max SH order 
%num_func:      Number of functions
%is_real:       Logical, if true, C is real, otherwise, C is complex
%is_real_dec:   Logical, if true, Y * C is real (Decoded random function is real in time-domain)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C              [(max_odr + 1)^2 x num_func] SH coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Mercator plot of random functions

% rng(234);
% max_odr = 3;
% num_func = 1;
% mode = 'mercator';
% C_cpx = sh_rand(max_odr, num_func, false, false);
% sh_plt(C_cpx, mode, false, 'title_name', 'Complex, Decode Complex');
% C_re = sh_rand(max_odr, num_func, true, false);
% sh_plt(C_re, mode, true, 'title_name', 'Real, Decode Complex');
% C_re_dec = sh_rand(max_odr, num_func, false, true);
% sh_plt(C_re_dec, mode, false, 'title_name', 'Complex, Decode Real');
% C_re_re_dec = sh_rand(max_odr, num_func, true, true);
% sh_plt(C_re_re_dec, mode, true, 'title_name', 'Real, Decode Real');

arguments
    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;
    num_func (1,1) double {mustBeNonnegative, mustBeInteger} = 1;
    is_real (1,1) logical = false;
    is_real_dec (1,1) logical = false;
end

if is_real
    C = randn((max_odr + 1)^2, num_func);
else
    C = randn((max_odr + 1)^2, num_func) + randn((max_odr + 1)^2, num_func) * 1i;
end

if is_real_dec
    C = (C + sh_conj(C)) / 2;
end