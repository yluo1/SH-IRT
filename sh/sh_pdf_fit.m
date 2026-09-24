function [C_pdf, err] = sh_pdf_fit(X, theta, phi, max_odr, is_real, mode, options)
%Fit spherical harmonic expansion probability density function to target
%evaluations X at spherical coordinates theta, phi

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X:             [N x M] N PDF evaluations of M functions, non-negative

%theta:         [N x 1] Co-latitude [0, pi]
%phi:           [N x 1] Azimuth [0, 2 * pi)

%max_odr:       Max SH order 

%is_real:       Logical, if true, evaluate real SH

%mode:          String, fitting method {'SqProjNNLS', 'SqProjQP', 'SqMagMS', 'SqMagSOMS'}
%                       'SqProjNNLS':       Solve for non-negative least squares weights 
%                                           of magnitude squared projection kernels, normalize to unity integral constraints in post
%                       'SqProjQP':         Solve for non-negative least squares weights, subject to unity integral constraints
%                                           of magnitude squared projection kernels 
%                       'SqMagMS':          Solve for magnitude square least squares, subject to unity integral constraints
%                                           with generic fmincon solver
%                       'SqMagSOMS':        Solve for sum-of-magnitude square least squares, subject to unity integral constraints
%                                           with semi-definite programming (cvx)

%options:       struct

%options.SqProjQP_options_quadprog:     struct, SqProjQP options for quadprog
%options.SqMagFmincon_C0:               [(floor(max_odr/2)+1)^2 x K]  SqMagFmincon initial guess for D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C_pdf:             [(max_odr + 1)^2 x M] SH coefficients
%err:               [1 x M] sum of square errors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Fit density to sample density

% max_odr = 8;
% N = (max_odr + 1)^2;
% is_real = false;
% rng(12463);
% dB_lim = [-60, 0];
% disp_phase = false;
% 
% [theta, phi] = sh_fib(N);
% %X = zeros(N, 1);  X(25) = 1; X(45) = 1;
% C_pdf_ref = sh_nrm(sh_msq(sh_rand(floor(max_odr/2), 1, is_real), is_real), 'Sum');
% X = real(sh_dec(C_pdf_ref, theta, phi, is_real));

%[C_pdf_SqProjNNLS, err_SqProjNNLS]         = sh_pdf_fit(X, theta, phi, max_odr, is_real, 'SqProjNNLS');
%[C_pdf_SqProjQP, err_SqProjQP]             = sh_pdf_fit(X, theta, phi, max_odr, is_real, 'SqProjQP');
%[C_pdf_SqMagFmincon, err_SqMagFmincon]     = sh_pdf_fit(X, theta, phi, max_odr, is_real, 'SqMagFmincon', 'SqProjQP_C0', sh_rand(floor(max_odr/2), 50, true));

%sh_plt(C_pdf_ref, 'mercator', is_real, 'dB_lim', dB_lim, 'disp_phase', disp_phase, 'title_name', 'PDF Ref.');

%sh_plt(C_pdf_SqProjNNLS, 'mercator', is_real, 'dB_lim', dB_lim, 'disp_phase', disp_phase, 'title_name', 'SqProjNNLS');
%err_SqProjNNLS

%sh_plt(C_pdf_SqProjQP, 'mercator', is_real, 'dB_lim', dB_lim, 'disp_phase', disp_phase, 'title_name', 'SqProjQP');
%err_SqProjQP

%sh_plt(C_pdf_SqMagFmincon, 'mercator', is_real, 'dB_lim', dB_lim, 'disp_phase', disp_phase, 'title_name', 'SqProjQP');
%err_SqMagFmincon


arguments
    X (:,:) double {mustBeNonnegative} = 0;

    theta (:,1) double = [0];
    phi   (:,1) double = [0];

    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;

    is_real (1,1) logical = false;

    mode (1,:) char {mustBeMember(mode, {'SqProjNNLS', 'SqProjQP', 'SqMagMS', 'SqMagSOMS'})} = 'SqProjQP';

    options.options_quadprog = optimoptions('quadprog');
    options.SqProjQP_C0 (:,:) double = [];
end

N = numel(theta);
assert(N == numel(phi), 'Size theta, phi mismatch');
assert(N == size(X, 1), 'Size theta, X mismatch');

M = size(X, 2);

%Check for duplicates
assert(size(unique([theta, phi], "rows"), 1) == N, '[theta, phi] must be unique');


C_pdf = zeros([(max_odr + 1)^2, M]);
err = zeros(1, M);

if strcmp(mode, 'SqProjNNLS')   %Non-negative least squares with post-normalization for unity integral

    C = sh_enc_proj_msq(max_odr, theta, phi, is_real); %[(max_odr + 1)^2 x N]
    A = real(sh_dec(C, theta, phi, is_real)); %[N x N] A(i,j) is contribution of kernel j to spherical coordinate i
    for m = 1:M
        [w_m, resnorm, residual] = lsqnonneg(A, X(:, m)); %[N x 1] non-negative weights per kernel expansion
        C_pdf(:, m) = sh_nrm(C * w_m, 'Sum');
    end

elseif strcmp(mode, 'SqProjQP') %Quadratic programming with linear constraints on unity integral

    C = sh_enc_proj_msq(max_odr, theta, phi, true); %[(max_odr + 1)^2 x N], real
    A = real(sh_dec(C, theta, phi, true)); %[N x N] A(i,j) is contribution of kernel j to spherical coordinate i
    for m = 1:M              
        [w_m, fval_m, exitflag_m] = quadprog((A')*A, - A*X(:, m), ...
            [], [], C(1,:), 1 / (2*sqrt(pi)), zeros(N,1), 1 ./ (C(1,:).' * 2 * sqrt(pi)), [], options.options_quadprog);

        C_pdf(:, m) = sh_nrm(C * w_m, 'Sum');
    end
    if ~is_real
        C_pdf = sh_re2cpx(C_pdf);
    end

elseif strcmp(mode, 'SqMagMS') %Magnitude square SH expansion least squares fitting 
                                    %with unit energy constraints on the squared function coefficient

    if isempty(options.SqProjQP_C0)
        C0 = sh_nrm(sh_enc_uni(floor(max_odr/2)), 'Sum'); %Uniform distribution
    else
        assert(size(options.SqProjQP_C0, 1) == (floor(max_odr/2) + 1)^2, 'Invalid options.SqProjQP_C0 size');
        C0 = options.SqProjQP_C0;
    end

    [C_pdf] = sh_fit_msq(X, theta, phi, max_odr, is_real, 'MS', ...
        'is_pdf', true, 'C0', C0);
    C_pdf = sh_nrm(C_pdf, 'Sum');
    
elseif strcmp(mode, 'SqMagSOMS') %Sum-of-magnitude square SH expansion least squares fitting

    [C_pdf] = sh_fit_msq(X, theta, phi, max_odr, is_real, 'SOMS', ...
        'is_pdf', true);
    C_pdf = sh_nrm(C_pdf, 'Sum');

else
    error('Unsupported mode');
end

%Compute errors
for m = 1:M
    err(m) = sum((real(sh_dec(C_pdf(:, m), theta, phi, is_real)) - X(:, m)).^2, 1);
end

