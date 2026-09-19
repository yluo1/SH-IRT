function [C, h_fig] = sh_enc_pist_sphere(max_odr, theta_center, phi_center, radius, theta_breadth, hz, r, options)
%Compute frequency response of a piston on a spherical baffle in the spherical harmonics domain

%Paper reference:
%Advanced Numerical Solutions. (n.d.). Exterior piston baffled in a sphere radiation problem (Coustyx Validation Manual). Ansol. ansol.us
%http://ansol.us/Products/Coustyx/Validation/Indirect/BoundaryExcitation/Sphere/SphericalBaffle/Downloads/dataset_description.pdf

%Implementation:    Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%max_odr:           Max SH order
%theta_center:      Spherical piston center's co-lattitude (radians)
%phi_center:        Spherical piston center's azimuth (radians)
%radius:            Spherical baffle's radius (meters)
%theta_breadth:     Spherical piston angular breadth (radians)

%hz:                [1 x N] Evaluation frequencies (Hz)
%r:                 Evaluation distance (Meters)

%options:                   Struct
%options.c:                 Speed of sound (meters / second)
%options.po:                Air density (kg / meter^3)
%options.u0:                Initial velocity (meters / second)
%options.normalize_onaxis:  Logical, if true, normalize to piston's center axis to unity

%options.enable_disp:       Logical, if true, plot frequency responses on horizontal plane

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:             [(max_odr+1)^2 x N] SH coefficients
%h_fig:         Figure handle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate sample frequency responses of exterior piston on sphere

%hz = logspace(log10(10), log10(24000), 512);

%C_ff = sh_enc_pist_sphere(40, pi/2, 0, 0.25, deg2rad(20), hz, 4, 'enable_disp', true);             % front firing
%C_sf = sh_enc_pist_sphere(40, pi/2, deg2rad(90), 0.25, deg2rad(20), hz, 4, 'enable_disp', true);   % side firing
%C_df = sh_enc_pist_sphere(40, pi, 0, 0.25, deg2rad(20), hz, 4, 'enable_disp', true);               % down firing

arguments

    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;
    theta_center (1,1) double = pi/2;
    phi_center (1,1) double = 0;
    radius (1,1) double {mustBePositive} = 1;
    theta_breadth (1,1) double {mustBePositive} = pi/4;
    
    hz (1,:) double {mustBeNonnegative} = 1000;
    r (1,1) double {mustBePositive} = 1;

    options.c (1,1) double {mustBePositive} = 343;
    options.po (1,1) double {mustBePositive} = 1.21;
    options.u0 (1,1) double {mustBePositive} = 1;
    options.normalize_piston_axis (1,1) logical = false;
    
    options.enable_disp (1,1) logical = false;
end

N = numel(hz);

k       = 2 * pi * hz / options.c;
kr      = k * r;
kradius = k * radius;
z_0     = options.po * options.c ;

b = complex(zeros(max_odr + 1, N));
for l = 0:max_odr

    if l == 0
        u_l = 0.5 * options.u0 * (1 - legendre_poly(l + 1, cos(theta_breadth)) );
    else
        u_l = 0.5 * options.u0 * (legendre_poly(l - 1, cos(theta_breadth)) - legendre_poly(l + 1, cos(theta_breadth)) );
    end
    A_l = -(1i * z_0) * (2 * l + 1) * u_l ./ ( l * (sh_hankel(l - 1, kradius, 'first') - sh_hankel(l + 1, kradius, 'first')) - sh_hankel(l + 1, kradius, 'first') );

    h = sh_hankel(l, kr, 'first');    
       
    b(l + 1, :) = (A_l .* h);
end
b = conj(b);

% Legenre addition theorem for centering spherical harmonics expansion
C = sh_val(max_odr, theta_center, phi_center, false)'; 

C = repmat(C, [1, N]);
for L = 0:max_odr
    idx = (L^2 + 1) : (L+1)^2;
    C(idx, :) = C(idx, :) .* b(L+1, :) * 4 * pi;
end
C(:, hz == 0) = 0;

% Normalize along piston axis direction
if options.normalize_piston_axis
    C(:, hz == 0) = 1;
    C = bsxfun(@rdivide, C, sh_dec(C, theta_center, phi_center, false));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_fig = [];
if options.enable_disp && coder.target('MATLAB')
    h_fig = sh_plt(C, 'horizontal', false, 'hz', hz);
end
