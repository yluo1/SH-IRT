function f = ft_exp_conv_opt_wrap(h, g)
%Codegen wrapper function ft_exp_conv_opt.m
%Convolution of FIR h by time-varying exponentiating FIR g

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%h:                 [1 x M] Fixed IR
%g:                 [1 x N] Exponentiating IR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%f:                 [1 x M * N]     Filtered IR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('ft_exp_conv_opt_wrap', '-o', 'ft/ft_exp_conv_opt_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare runtime and accuracy between Matlab and Mex

% rng(4123);
% Fs = 48000;
% T = 1;
% M = ceil(T * Fs);
% h = randn(1, M);
% g = ft_two_tap_FIR(0, -1); %Two tap
% tic; f_mat = ft_exp_conv_opt(h, g); duration_mat = toc
% tic; f_mex = ft_exp_conv_opt_mex(h, g); duration_mex = toc
% err = norm(f_mat - f_mex)

% g8 = ft_exp_conv_g(g, 7, 100);
% tic; f_mat = ft_exp_conv_opt(h, g8); duration_mat = toc
% tic; f_mex = ft_exp_conv_opt_mex(h, g8); duration_mex = toc
% err = norm(f_mat - f_mex)

arguments
    h (1,:) double = [1 0];
    g (1,:) double = 1;
end

f = ft_exp_conv_opt(h, g);
