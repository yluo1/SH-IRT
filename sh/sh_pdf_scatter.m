function D_pdf = sh_pdf_scatter(C_pdf, mode, lambda, is_real, options)
%Transport the spherical harmonic density function towards uniform distribution
%via attenuation of positive order SH coefficients

%mode:      'lin'  Linear interpolation between input and uniform distributions
%D_pdf(1, :)     = C_pdf(1, :);
%D_pdf(2:end, m) = lambda^(m-1) * C_pdf(2:end, m);

%mode:      'exp'   Exponential attenuation w.r.t. SH degree l
%idx = (2*l+1)^2 + m
%D_pdf(idx, n) = lambda^(options.exp_k * l * (n - 1)) * C_pdf(idx, n)

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
%C:                 [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%mode               String, attenuation method  {'lin', 'exp'}
%lambda:            Decay factor between [0, 1]
%is_real:           Logical, if true, evaluate real SH

%options:           struct
%options.exp_k:     Exponential decay factor per SH degree l, non-negative 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%D:                 [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Generate sequence of exponential decaying pdfs from Dirac expansion

%max_odr = 12;
%is_real = true;
%M = 8000;
%Fs = 16000;

%C_pdf = sh_nrm(sh_enc_rbf('SqExp', max_odr, pi/2, 0, 0.25, is_real), 'Sum');

%C_pdf = C_pdf * ones(1, M);
%D_lin = sh_pdf_scatter(C_pdf, 'lin', 0.999, is_real);
%D_exp_1 = sh_pdf_scatter(C_pdf, 'exp', 0.999, is_real, 'exp_k', 1);
%D_exp_half = sh_pdf_scatter(C_pdf, 'exp', 0.999, is_real, 'exp_k', 0.5);

%dB_lim = [-40, 20];
%t = (0:(M-1)) / Fs;
%sh_plt(C_pdf, 'horizontal', is_real, 'disp_phase', false, 'title_name', 'Original', 'dB_lim', dB_lim, 't', t);
%sh_plt(D_lin, 'horizontal', is_real, 'disp_phase', false, 'title_name', 'Linear Scattering', 'dB_lim', dB_lim, 't', t);
%sh_plt(D_exp_1, 'horizontal', is_real, 'disp_phase', false, 'title_name', 'Exponential k = 1 Scattering', 'dB_lim', dB_lim, 't', t);
%sh_plt(D_exp_half, 'horizontal', is_real, 'disp_phase', false, 'title_name', 'Exponential k = 0.5 Scattering', 'dB_lim', dB_lim, 't', t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Generate sequence of exponential decaying pdfs from DiracRandom expansion

%max_odr = 12;
%is_real = true;
%M = 8000;
%Fs = 16000;

%rng(123);
%C_pdf = sh_pdf_preset('FwdCenter', 'DiracRandom', max_odr, is_real);
%C_pdf = C_pdf * ones(1, M);
%D_lin = sh_pdf_scatter(C_pdf, 'lin', 0.999, is_real);
%D_exp_1 = sh_pdf_scatter(C_pdf, 'exp', 0.999, is_real, 'exp_k', 1);
%D_exp_half = sh_pdf_scatter(C_pdf, 'exp', 0.999, is_real, 'exp_k', 0.5);

%dB_lim = [-24, -6];
%t = (0:(M-1)) / Fs;
%sh_plt(C_pdf, 'horizontal', is_real, 'disp_phase', false, 'title_name', 'Original', 'dB_lim', dB_lim, 't', t);
%sh_plt(D_lin, 'horizontal', is_real, 'disp_phase', false, 'title_name', 'Linear Scattering', 'dB_lim', dB_lim, 't', t);
%sh_plt(D_exp_1, 'horizontal', is_real, 'disp_phase', false, 'title_name', 'Exponential k = 1 Scattering', 'dB_lim', dB_lim, 't', t);
%sh_plt(D_exp_half, 'horizontal', is_real, 'disp_phase', false, 'title_name', 'Exponential k = 0.5 Scattering', 'dB_lim', dB_lim, 't', t);

arguments
     C_pdf (:,:) double {coder.mustBeComplex} = complex(0);
     mode (1,:) char {mustBeMember(mode, {'lin', 'exp'})} = 'exp'
     lambda (1,1) double {mustBeNonnegative} = 0.99;
     is_real (1,1) logical = false;

     options.exp_k (1,1) double {mustBeNonnegative} = 0.5;
end

assert(sh_pdf_check(C_pdf, is_real), 'C_pdf is invalid density');

[N_C, M] = size(C_pdf);
P = sqrt(N_C) - 1;
assert(P - floor(P) == 0, 'C_pdf invalid size');

assert(lambda <= 1 && lambda >= 0, 'lambda must be between [0, 1]')

D_pdf = C_pdf;

if strcmp(mode, 'lin')

    D_pdf(2:end, :) = bsxfun(@times, D_pdf(2:end, :), lambda.^(0:(M-1)) );

elseif strcmp(mode, 'exp')

    for l = 0:P
        idx = (l^2 + 1) : (l+1)^2;
        D_pdf(idx, :) = bsxfun(@times, D_pdf(idx, :), lambda.^(options.exp_k * l * (0:(M-1))  ) );
    end
    
else
    error('Unsupported mode')
end
