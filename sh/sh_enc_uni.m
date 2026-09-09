function C = sh_enc_uni(max_odr, enable_disp)
%Spherical harmonic expansion of unity function
%f(\Omega) = 1 \forall \Omega

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%max_odr:       Max SH order 
%enable_disp:   Logical, if true, plot function expansion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:             [(max_odr + 1)^2 x N] SH coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_enc_uni', '-o', 'sh/sh_enc_uni_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Plot and decode unity function

% max_odr = 3;
% rng(343);
% N_dec = 100;
% C = sh_enc_uni(max_odr, true);
% f = sh_dec(C, rand(N_dec,1) * pi, rand(N_dec,1));
% err = norm(f - 1)

arguments
    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;
    enable_disp (1,1) logical = false;
end

C = complex(zeros([(max_odr + 1)^2, 1]));
C(1) = 2 * sqrt(pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if enable_disp && coder.target('MATLAB')

    %Mercator plot
    sh_plt(C, 'mercator', false);

end