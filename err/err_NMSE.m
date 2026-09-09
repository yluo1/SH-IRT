function err = err_NMSE(obs, ref)
%Normalized mean squared error
%sum((obs - ref).^2) / sum(ref.^2)

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%obs:         [N x 1] or [1 x N] Observed 
%ref:         [N x 1] or [1 x N] Reference

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%err:       Error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(numel(obs) == numel(ref), 'a, b size mismatch');

err = sum((obs(:) - ref(:)).^2) / sum(ref(:).^2);
