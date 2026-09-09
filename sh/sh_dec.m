function f = sh_dec(C, theta, phi, is_real, mode)
%Evaluate (decode) spherical harmonic expansion at select spherical coordinates

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input 
%C:             [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)

%theta:         [N x 1]  Co-latitude [0, pi]
%phi:           [N x 1]  Azimuth [0, 2 * pi)

%is_real:       Logical, if true, evaluate real SH

%mode:          String, decoding method {'eval', 'pinv'}
%               'eval':     Evaluate spherical harmonic expansion Y(theta, phi) * C 
%               'pinv':     Pseudo-inverse (least-squares equivalent) Y(theta, phi) \ C

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%f:             [N x M] SH decoded at N points across M functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_dec', '-o', 'sh/sh_dec_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare Matlab versus Mex

% rng(6451);
% P = 5;
% M = 6;
% C = sh_rand(P, M);

% N = 200;
% theta = rand(N, 1) * pi;
% phi   = rand(N, 1) * 2 * pi;
% mode = 'eval';

% tic; f_cpx_mat = sh_dec(C, theta, phi, false, mode); toc_cpx_mat = toc
% tic; f_cpx_mex = sh_dec_mex(C, theta, phi, false, mode); toc_cpx_mex = toc
% err_cpx = norm(f_cpx_mat - f_cpx_mex)

% tic; f_real_mat = sh_dec(real(C), theta, phi, true, mode); toc_real_mat = toc
% tic; f_real_mex = sh_dec_mex(real(C), theta, phi, true, mode); toc_real_mex = toc
% err_real = norm(f_real_mat - f_real_mex)

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
    
    theta (:,1) double = [0];
    phi   (:,1) double = [0];

    is_real (1,1) logical = false;

    mode (1,:) char {mustBeMember(mode, {'eval', 'pinv'})} = 'eval';
end

[P, M] = size(C);
P = sqrt(P) - 1;

assert(P - floor(P) == 0, 'Invalid size C');
assert(numel(theta) == numel(phi), 'theta, phi size mismatch');
N = numel(theta);

if strcmp(mode, 'eval')
    %Evaluate at spherical coordinates
    Y = sh_val(P, theta, phi, is_real);  %[N x (P + 1)^2] 
    f = Y * C;

elseif strcmp(mode, 'pinv')
    %Least-squares (Pseudo-inverse equivalent)
    assert(N >= (P+1)^2, 'Requires N >= (P+1)^2');
    Y = sh_val(P, theta, phi, is_real);  %[N x (P + 1)^2]
    f = Y \ C;
    
else
    error('Unsupported mode');
end