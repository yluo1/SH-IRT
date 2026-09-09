function C = sh_rand_proj(max_odr, num_func, is_real, is_mag_sq)
%Generate random projection SH expansions 

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%max_odr:       Max SH order 
%num_func:      Number of functions
%is_real:       Logical, if true, C is real, otherwise, C is complex
%is_mag_sq:     Logical, if true, generate magnitude squared projections at floor(max_odr/2) order

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C              [(max_odr + 1)^2 x num_func] SH coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Mercator plot of sum of random projections

%rng(2314);
%max_odr = 4;
%is_real = false;
%C_1 =  sh_nrm(sum(sh_rand_proj(max_odr, 1, is_real), 2), 'Sum');
%C_10 = sh_nrm(sum(sh_rand_proj(max_odr, 10, is_real), 2), 'Sum');
%C_10_msq = sh_nrm(sum(sh_rand_proj(max_odr, 10, is_real, true), 2), 'Sum');

%sh_plt(C_1, 'mercator', is_real, 'dB_lim', [-60, 0]);
%sh_plt(C_10, 'mercator', is_real, 'dB_lim', [-60, 0]);
%sh_plt(C_10_msq, 'mercator', is_real, 'dB_lim', [-60, 0]);

arguments
    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;
    num_func (1,1) double {mustBeNonnegative, mustBeInteger} = 1;
    is_real (1,1) logical = false;
    is_mag_sq (1,1) logical = false;
end

if is_mag_sq
    C = sh_msq(sh_enc_proj(floor(max_odr/2), acos(2 * rand(num_func, 1) - 1) , rand(num_func, 1) * 2 * pi, is_real), is_real);
    C = [C; zeros((max_odr + 1)^2 - size(C, 1), num_func)];
else
    C = sh_enc_proj(max_odr, acos(2 * rand(num_func, 1) - 1) , rand(num_func, 1) * 2 * pi, is_real);
end