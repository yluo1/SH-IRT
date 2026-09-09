function D = sh_nrm(C, mode)
%N3D Spherical harmonic expansion normalizations

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:             [(P + 1)^2 x M]  SH coefficients (max order P of M number of functions)

%mode:          String, normalization method
%               'Sum'            Summation (Integral of f(\Omega) d\Omega over sphere = 1), required for probability density functions
%               'Pow'            Total power (Integral of |f(\Omega)|^2 d\Omega over sphere = 1)
%               'PowAvg'         Power average (Integral of |f(\Omega)|^2 d\Omega / (4 * pi) over sphere = 1)

%               'N3D'            Full normalization (identity)
%               'SN3D'           Schmidt semi-normalization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ouput
%D:             [(P + 1)^2 x M] SH coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_nrm', '-o', 'sh/sh_nrm_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare normalization methods of random function

% rng(2351);
% max_odr = 5;
% C = sh_rand(max_odr);
% C_sum = sh_nrm(C, 'Sum');
% sh_int(C_sum, 'Sum')
% C_pow = sh_nrm(C, 'Pow');
% sh_int(C_pow, 'SumSq')
% C_pow_avg = sh_nrm(C, 'PowAvg');
% sh_int(C_pow_avg, 'PowAvg')
% varargin = {'mercator', true, 'dB_lim', [-20, 20]};
% sh_plt(C, varargin{:}, 'title_name', 'Original');
% sh_plt(C_sum, varargin{:}, 'title_name', 'Normalized Sum');
% sh_plt(C_pow, varargin{:}, 'title_name', 'Normalized Total Power');
% sh_plt(C_pow_avg, varargin{:}, 'title_name', 'Normalized Power Average');

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
    mode (1,:) char {mustBeMember(mode, {'Sum', 'Pow', 'PowAvg', 'N3D', 'SN3D'})} = 'Sum';
end

[P, M] = size(C);
P = sqrt(P) - 1;
assert(P - floor(P) == 0, 'Invalid size C');

%Normalization methods
if strcmp(mode, 'Sum') %Summation

    D = bsxfun(@rdivide, C, sh_int(C, 'Sum') );

elseif strcmp(mode, 'Pow') %Total power 

    D = bsxfun(@rdivide, C, sqrt(sh_int(C, 'SumSq')) );    

elseif strcmp(mode, 'PowAvg') %Power average

    D = bsxfun(@rdivide, C, sqrt(sh_int(C, 'PowAvg')) );    

elseif strcmp(mode, 'SN3D')  %Schmidt semi-normalization

    D = C;
    for l = 0:P
        idx = (l^2 + 1) : (l+1)^2;
        D(idx, :) = C(idx, :) ./ sqrt(2 * l + 1);              
    end

elseif strcmp(mode, 'N3D')    %Full normalization (identity)

    D = C;

else

    D = C;
    error('Unsupported mode');

end