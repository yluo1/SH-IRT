function [f, g] = sh_msq_err(C, X, Y)
%Magnitude squared SH expansion error
%f(C) = abs(Y(theta, phi) * C)^2 - X(theta, phi)

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C: [(P + 1)^2 x 1]             Spherical harmonic expansion coefficients
%X: [N x 1]                     N target observations, non-negative
%Y: [N x (P + 1)^2]             Spherical harmonic bases evaluations 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%f:     Scalar, sum of square errors
%g:     [(P + 1)^2 x 1] gradient g(C) = df(C)/dC

arguments
     C (:,1) double = 0;     
     X (:,1) double = 1;
     Y (:,:) double = 1;
end

f_c = Y * C;                     %[N x 1]
fc_c_mag_sq = f_c .* conj(f_c);  %[N x 1]
err = fc_c_mag_sq - X;           %[N x 1]
f = sum(err.^2);

%f = (|Y*C|.^2 - X)' * (|Y*C|.^2 - X)
%  = err(C)' * err(C)

%derr(C)/dC = d(Y*C)/dC .* (Y*C) + (Y*C) .* d(Y*C)/dC,  %Product rule, Hadamard product commutes
%           = 2 * (Y*C) .* d(Y*C)/dC
%           = 2 * diag(Y*C) * Y

%df/dC = 2 * derr/dC' * err 
%      = 4 * (diag(Y*C) * Y)' * err
%      = 4 * Y' * diag(Y*C) * err
%      = 4 * Y' * f_c .* err

g = 4 * Y' * (f_c .* err);

