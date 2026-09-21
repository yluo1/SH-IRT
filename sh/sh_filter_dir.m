function E = sh_filter_dir(C, D, is_real, mode)
% Directional filtering via SH multiplication

%f(theta, phi) = C.' * Y(theta, phi)
%g(theta, phi) = D.' * Y(theta, phi)

%f(theta, phi) * g(theta, phi) = E.' * Y(theta, phi)

%Express output SH coefficients in terms of transfer matrix A vector product:
%E(:, m) = A(:, :, m) * D(:, m) forall 1 <= m <= M

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:             [(P_C + 1)^2 x M]   SH coefficients (max order P_C of M number of functions)
%D:             [(P_D + 1)^2 x M]   SH coefficients (max order P_D of M number of functions)

%is_real:       Logical, if true, evaluate real SH

%mode:          String, muliplication method {'CGC', 'SHT'}
%                   'CGC':   Clebsch-Gordan coefficients
%                   'SHT':   Spherical harmonic transform (inverse -> prod -> forward)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ouput
%E:             [(P_C + P_D + 1)^2 x M]                     SH coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Compare CGC and SHT modes

% rng(2134);
% P_C = 7;
% P_D = 7;
% 
% theta = pi/2;
% phi = 0;
% ell = 0.5;
% is_real = false;
% D = sh_enc_rbf('SqExp', P_D, theta, phi, ell, is_real, false);
% C = sh_rand(P_C, 1);
% tic; E_CGC = sh_filter_dir(D, C, is_real, 'CGC');  duration_CGC = toc
% tic; E_SHT = sh_filter_dir(D, C, is_real, 'SHT');  duration_SHT = toc

% varargin = {'dB_lim', [-32, 18], 'disp_phase', false,  'fontsize', 20};

% h_D = sh_plt(D, 'mercator', false, 'title_name', 'Directional Function', varargin{:}, 'disp_cb', false, 'disp_ylabel', true, 'disp_yticks', true);
% h_C = sh_plt(C, 'mercator', false, 'title_name', 'Random Function', varargin{:}, 'disp_cb', false, 'disp_ylabel', false, 'disp_yticks', false);
% h_E_CGC = sh_plt(E_CGC, 'mercator', false, 'title_name', 'Product Function (CGC)', varargin{:}, 'disp_cb', true, 'disp_ylabel', false, 'disp_yticks', false);
% h_E_SHT = sh_plt(E_SHT, 'mercator', false, 'title_name', 'Product Function (SHT)', varargin{:}, 'disp_cb', true, 'disp_ylabel', false, 'disp_yticks', false);

% err = norm(E_CGC - E_SHT)

arguments
    C (:,:) double {coder.mustBeComplex} = complex(0)
    D (:,:) double {coder.mustBeComplex} = complex(0)

    is_real (1,1) logical = false;
    
    mode (1,:) char {mustBeMember(mode, {'CGC', 'SHT'})} = 'SHT';
end

[P_C, M_C] = size(C);
[P_D, M_D] = size(D);
P_C = sqrt(P_C) - 1;
P_D = sqrt(P_D) - 1;
P_E = P_C + P_D;

assert(M_C == M_D, 'M not matched for C, D');
assert(P_C - floor(P_C) == 0, 'Invalid size C');
assert(P_D - floor(P_D) == 0, 'Invalid size D');

N_E = (P_E + 1)^2;

if strcmp(mode, 'CGC')

    if is_real % Convert to complex form
        C = sh_re2cpx(C);
        D = sh_re2cpx(D);
        E_cpx = sh_mul(C, D);
        E = sh_cpx2re(E_cpx); %Convert back to real form
    else
        E = sh_mul(C, D);
    end

elseif strcmp(mode, 'SHT')

    [theta, phi] = sh_fib(N_E);
    
    f = sh_dec(C, theta, phi, is_real) .* sh_dec(D, theta, phi, is_real);

    Y = sh_val(P_E, theta, phi, is_real);
    E = Y \ f;
    
else
    error('Unknown mode')
end

