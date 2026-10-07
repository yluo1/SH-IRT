function err = err_SNMSE(a, b)
%Symmetric normalized mean squared error
%sum(abs(a(:) - b(:)).^2) / ( (sum(abs(a(:)).^2) + sum(abs(b(:)).^2)) / 2 )

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%a:         [N x 1] or [1 x N] 
%b:         [N x 1] or [1 x N]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%err:       Error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(numel(a) == numel(b), 'a, b size mismatch');

err = sum(abs(a(:) - b(:)).^2) / ( (sum(abs(a(:)).^2) + sum(abs(b(:)).^2)) / 2 );