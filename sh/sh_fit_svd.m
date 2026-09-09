function [C, err, trunc_percent] = sh_fit_svd(X, theta, phi, max_odr, is_real, trunc_frac)
%Least-squares spherical harmonic basis fit min ||Y*C - X||^2 
%via truncated singular value decomposition:

%singular_values(singular_values <= max(singular_values) * trunc_frac) = 0

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X:             [N x M] N measurements of M functions

%theta:         [N x 1] Co-latitude [0, pi]
%phi:           [N x 1] Azimuth [0, 2 * pi)

%max_odr:       Max SH order 

%is_real:       Logical, if true, evaluate real SH

%trunc_frac:    Fraction of largest singular values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:                 [(max_odr + 1)^2 x M] SH coefficients
%err:               Scalar, norm(Y*C - X);
%trunc_percent:     Percentage of singular values truncated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_fit_svd', '-o', 'sh/sh_fit_svd_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Least-squares SH fit to random function distributed along 
%spherical Fibonnaci points with varying truncation fractions

% rng(441);
% max_odr = 10;
% is_real = false;
% N_pts = (max_odr + 1)^2;
% X = randn(N_pts, 1) + randn(N_pts, 1) * 1i;
% [theta, phi] = sh_fib(N_pts);

% [C_trunc_0, err_trunc_0, trunc_percent_0] = sh_fit_svd(X, theta, phi, max_odr, is_real, 0); err_trunc_0
% sh_plt(C_trunc_0); trunc_percent_0

% [C_trunc_01, err_trunc_01, trunc_percent_01] = sh_fit_svd(X, theta, phi, max_odr, is_real, 0.01); err_trunc_01
% sh_plt(C_trunc_01); trunc_percent_01

% [C_trunc_50, err_trunc_50, trunc_percent_50] = sh_fit_svd(X, theta, phi, max_odr, is_real, 0.50); err_trunc_50
% sh_plt(C_trunc_50); trunc_percent_50

% [C_trunc_90, err_trunc_90, trunc_percent_90] = sh_fit_svd(X, theta, phi, max_odr, is_real, 0.90); err_trunc_90
% sh_plt(C_trunc_90); trunc_percent_90

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Least-squares SH fit to radial basis function along 
%spherical Fibonnaci points with varying max_odr

% rng(441);
% max_odr = 10;
% is_real = false;
% N_pts = (max_odr + 1)^2;
% X = randn(N_pts, 1) + randn(N_pts, 1) * 1i;
% [theta, phi] = sh_fib(N_pts);

% [C_max_odr_9, err_max_odr_9] = sh_fit_svd(X, theta, phi, 9, is_real, 0); 
% sh_plt(C_max_odr_9); err_max_odr_9

% sh_plt(sh_val(max_odr, theta, phi, is_real) \ X);

% [C_max_odr_11, err_max_odr_11] = sh_fit_svd(X, theta, phi, 11, is_real, 0); 
% sh_plt(C_max_odr_11); err_max_odr_11


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare Matlab versus Mex

% rng(441);
% max_odr = 10;
% is_real = false;
% N_pts = (max_odr + 1)^2;
% X = randn(N_pts, 1) + randn(N_pts, 1) * 1i;
% [theta, phi] = sh_fib(N_pts);

% C_trunc_mat = sh_fit_svd(X, theta, phi, max_odr, is_real, 0.90); 
% C_trunc_mex = sh_fit_svd_mex(X, theta, phi, max_odr, is_real, 0.90); 
% err = norm(C_trunc_mat - C_trunc_mex)

arguments
    X (:,:) double {coder.mustBeComplex} = complex(0);

    theta (:,1) double = [0];
    phi   (:,1) double = [0];

    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;

    is_real (1,1) logical = false;

    trunc_frac (1,1) double {mustBeNonnegative, mustBeLessThanOrEqual(trunc_frac, 1)} = 0;
end

%Compute SH bases
Y = sh_val(max_odr, theta, phi, is_real);  %[N x (max_odr + 1)^2]

%Compute singular value decomposition
[U, S, V] = svd(Y);
s_list = diag(S);   %List of singular values
N_s_list = numel(s_list);

%Compute pseudo-inverse singular values
inv_s_list = 1 ./ s_list;

%Truncate singular values <= frac_trunc * largest singular value 
trunc_mask = s_list <= (max(s_list) * trunc_frac);
inv_s_list(trunc_mask) = 0;
trunc_percent = sum(trunc_mask) / N_s_list * 100;

inv_S = S;
inv_S(1:N_s_list, 1:N_s_list) = diag(inv_s_list);
C = V * inv_S' * U' * X;

%Compute error
err = norm(Y * C - X);