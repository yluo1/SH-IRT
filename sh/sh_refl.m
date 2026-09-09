function D = sh_refl(C, mode, is_real)
%Reflect spherical harmonic expansion over a cardinal plane

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:        [(P+1)^2 x M] SH coefficients (max order P of M number of functions)
%mode:     String, plane type {'xz', 'median', 'sagittal', 'xy', 'transverse', 'yz', 'frontal'}
%is_real:  Logical, if true, evaluate real SH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%D:         [(P+1)^2 x M] SH coefficients, P max-order, M functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Reflect sample SqExp SH expansion over cardnal planes

%max_odr = 3;
%theta = deg2rad(45);
%phi = deg2rad(45);
%ell = 1;
%is_real = false;
%C = sh_enc_rbf('SqExp', max_odr, theta, phi, ell, is_real);
%D_med = sh_refl(C, 'median', is_real);
%D_trans = sh_refl(C, 'transverse', is_real);
%D_front = sh_refl(C, 'frontal', is_real);

%sh_plt(C, 'mercator', is_real, 'title_name', 'Original');
%sh_plt(D_med, 'mercator', is_real, 'title_name', 'Reflected over Median');
%sh_plt(D_trans, 'mercator', is_real, 'title_name', 'Reflected over Transverse Plane');
%sh_plt(D_front, 'mercator', is_real, 'title_name', 'Reflected over Frontal Plane');

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
    mode (1,:) char {mustBeMember(mode, {'xz', 'median', 'sagittal', 'xy', 'transverse', 'yz', 'frontal'})} = 'median';
    is_real (1,1) logical = false;  
end

if strcmp(mode, 'xz') || strcmp(mode, 'median')  || strcmp(mode, 'sagittal')

    D = conj(sh_conj(C));

elseif strcmp(mode, 'xy')  || strcmp(mode, 'transverse')

    D = sh_rot(conj(sh_conj(sh_rot(C, [0 0 pi/2], is_real))), [0 0 -pi/2], is_real);
        
elseif strcmp(mode, 'yz')  || strcmp(mode, 'frontal')

    D = sh_rot(conj(sh_conj(sh_rot(C, [pi/2 0 0], is_real))), [-pi/2 0 0], is_real);
    
else
    
    error('Unsupported plane');

end