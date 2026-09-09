function [C, t] = sh_ism_wrap(max_src_odr, max_rec_odr, is_real, T, s, r, l, gamma_pos, gamma_neg, ...
    Fs, kernel_sample_width, jitter_coord_bnd, jitter_srand)
%Codegen wrapper function around sh_ism.m

%Generate all pair-wise combinations of spherical harmonic expansions
%between acoustic source and receiver in image-source-receiver room model.
%See paper reference:

%Luo, Y. and Kim, W., 2021, January. 
%Fast source-room-receiver acoustics modeling.
%In 2020 28th European Signal Processing Conference (EUSIPCO) (pp. 51-55). IEEE.

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%max_src_odr:   Max SH order of source expansions
%max_rec_odr:   Max SH order of receiver expansions
%is_real:       Logical, if true, C is real, otherwise, C is complex, non-negative

%T:             Scalar, time (sec)
%s:             [1 x 3] Source offset from origin (meters)
%r:             [1 x 3] Receiver offset from origin (Meters)
%l:             [1 x 3] Orthotope dimensions (meters)
%gamma_pos:     [L x 3] Filter (L-tap length) reflection IR for wall on +axis
%gamma_neg:     [L x 3] Filter (L-tap length) reflection IR for wall on -axis

%Fs:                    Sampling rate
%kernel_sample_width:   Lanczos kernel half-window size
%jitter_coord_bnd:      [1 x 2] (min, max) bounds on jitter to image-source and receiver coordinates
%jitter_srand:          Random seed for jitter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:     [(max_src_odr + 1)^2 x (max_rec_odr + 1)^2 x M]
%t:     [1 x M] time (sec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_ism_wrap', '-o', 'sh/sh_ism_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare Matlab with Mex accuracy and runtime

% max_src_odr = 1;
% max_rec_odr = 5;
% is_real = true;
% T = 0.2;
% s = [1 2 1];
% r = [2 0.7 1];
% l = [5 6 3];    
% gamma_pos = [0.8, 0.7, 0.5; 0.2, 0.1, 0.3];
% gamma_neg = [0.9, 0.6, 0.5; 0.1, 0.1, 0.2];

%Fs = 48000;
%kernel_sample_width = 10;
%jitter_coord_bnd = [1 1] * 1e-3;
%jitter_coord_bnd = [0, 0];
%jitter_srand = 6452;

%tic;
% C_mat = sh_ism(max_src_odr, max_rec_odr, is_real, T, s, r, l, gamma_pos, gamma_neg, ...
%               'Fs', Fs, 'kernel_sample_width', kernel_sample_width, 'jitter_coord_bnd', jitter_coord_bnd, 'jitter_srand', jitter_srand);
% duration_mat = toc

%tic;
% C_mex = sh_ism_mex(max_src_odr, max_rec_odr, is_real, T, s, r, l, gamma_pos, gamma_neg, ...
%                   Fs, kernel_sample_width,  jitter_coord_bnd,  jitter_srand);
% duration_mex = toc

%err = norm(C_mat(:) - C_mex(:))

arguments
    max_src_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 0;
    max_rec_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 3;
    is_real (1,1) logical = false;

    T (1,1) double {mustBeNonnegative} = 0.2;
    s (1,3) double = [1 2 1];
    r (1,3) double = [2 1 1];
    l (1,3) double {mustBePositive} = [5 6 3];
    
    gamma_pos (:,3) double = [0.8, 0.8, 0.8];
    gamma_neg (:,3) double = [0.8, 0.8, 0.8];


    Fs (1,1) double {mustBePositive} = 48000;
    kernel_sample_width  (1,1) double {mustBePositive, mustBeInteger}  = 10;
    jitter_coord_bnd  (1,2) double {mustBeNonnegative} = [0, 0];
    jitter_srand (1,1) double {mustBeInteger} = 6452;

end


[C, t] = sh_ism(max_src_odr, max_rec_odr, is_real, T, s, r, l, gamma_pos, gamma_neg, ...
    'Fs', Fs, 'kernel_sample_width', kernel_sample_width, 'jitter_coord_bnd', jitter_coord_bnd, 'jitter_srand', jitter_srand);

