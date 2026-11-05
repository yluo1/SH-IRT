function K = sh_cov(C, D_pdf, is_real, options)
% Compute spatial covariance between spherical harmonic expansions in C over density D_pdf in the far-field

% SH expansion: f_C(theta, phi) = Y(theta, phi) * C    row-vector [1 x N]
% K = E_{D_pdf}[f_C.' * conj(f_C)]
% K_{ij} = \int_{\Omega} Y(\Omega) * C(:, i) * conj( Y(\Omega) * C(:, j) ) * Y(\Omega) * D_pdf d\Omega for uniform grid \Omega

% mode = 'CGC':     Compute analytic covariance via Clebsch-Gordan coefficients (SH product operator) for common centers-of-expansions
% mode = 'SHT':     Compute numeric covariance via quadrature over uniform spherical coordinate grid.
%                   Supports non-common centers-of-expansion direction and delay compensation (require both options.SHT_centers and options.SHT_wavnum)

% Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:             [(L_C + 1)^2 x N x M] SH coefficients (max order L_C of N number of expansions, M functions)

%D_pdf:         [(L_D + 1)^2 x 1] or [(L_D + 1)^2 x M]       Density SH coefficients (max order L_D) across all/each of M functions

%is_real:       Logical, if true, evaluate real SH

%options.mode:          String, integration (SH multiplication) method {'CGC', 'SHT'}
%                           'CGC':   Clebsch-Gordan coefficients (SH product operator)
%                           'SHT':   Spherical harmonic transform (inverse SHT -> numeric integration)

%options.center_mean:       Logical, if true, subtract mean from covariance matrix
%options.SHT_oversample:    Oversample number of spherical coordinate direction by this factor for mode = 'SHT'
%options.SHT_maxsample:     Max number of spherical coordinate direction for mode = 'SHT'

%options.SHT_center:        [] or [N x 3]   Non-common centers-of-expansion coordinates [x, y, z]
%options.SHT_wavnum:        [] or [1 x M]   Angular wavenumber kappa (e.g. kappa = 2 * pi * hz / c for propagation velocity c)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%K:             [N x N x M]     Covariance matrix [N x N] per function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:  Compare CGC to SHT methods

% is_real = false;
% max_odr = 6;
% D_pdf = sh_pdf_preset('FwdCenter', 'SqExp', max_odr, is_real, 'enable_disp', true);
% sh_plt(D_pdf, 'mercator', is_real, 'disp_phase', false, 'title_name', 'Density');

%C_left  = sh_enc_pist_sphere(max_odr, pi/2, deg2rad(10), 0.25, deg2rad(20), 400, 4, 'enable_disp', true);
%C_right = sh_enc_pist_sphere(max_odr, pi/2, deg2rad(-10), 0.25, deg2rad(20), 400, 4, 'enable_disp', true);
%C = [C_left, C_right];

%center_mean = true;
%K_CGC = sh_cov(C, D_pdf, is_real, 'mode', 'CGC', 'center_mean', center_mean)
%K_SHT = sh_cov(C, D_pdf, is_real, 'mode', 'SHT', 'center_mean', center_mean)
%err = err_NMSE(K_SHT, K_CGC)

% C_left   = sh_enc_uni(max_odr);
% C_right  = sh_enc_uni(max_odr) .* 1i;
% C = [C_left, C_right];

%center_mean = false;
%K_CGC = sh_cov(C, D_pdf, is_real, 'mode', 'CGC', 'center_mean', center_mean)
%K_SHT = sh_cov(C, D_pdf, is_real, 'mode', 'SHT', 'center_mean', center_mean)
%err = err_NMSE(K_SHT, K_CGC)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare SHT for centered versus non-centered expansions

% Fs = 8000;
% is_real = false;
% max_odr = 12;
% N_taps = 128;
% w = linspace(0, 2 * pi, N_taps + 1); w = w(1:end-1);
% t = (0:(N_taps-1))/Fs;
% omega = w * Fs;
% freq = omega / (2 * pi);
% v = [- 0.1/1000 * 343, 0, 0];

% D_pdf = sh_pdf_preset('FwdCenter', 'SqExp', max_odr, is_real);
% C = sh_enc_pist_sphere(max_odr, pi/2, 0, 0.25, deg2rad(20), freq, 4, 'normalize_piston_axis', true);
% C_translate = sh_translate(C, omega, is_real, v, 'enable_disp', true, 'disp_dB_lim', [-20, 6]); % 1/2 ms in -x direction
% C_resize = sh_resize(C, sqrt(size(C_translate, 1)) - 1 );
% K_translate = sh_cov(permute(cat(3, C_resize, C_translate), [1, 3, 2]), D_pdf, is_real);
% K_noncenter = sh_cov(permute(cat(3, C_resize, C_resize), [1, 3, 2]), D_pdf, is_real, 'SHT_center', [0 0 0; v], 'SHT_wavnum', omega / 343);
% err = err_NMSE(K_translate, K_noncenter)

arguments
    C (:,:,:) double {coder.mustBeComplex} = complex(0);
    D_pdf (:,:) double {coder.mustBeComplex} = sh_pdf_preset('FwdCenter', 'Uniform', 0);
    is_real (1,1) logical = false;
    
    options.mode (1,:) char {mustBeMember(options.mode, {'CGC', 'SHT'})} = 'SHT';
    options.center_mean (1,1) logical = false;

    options.SHT_oversample (1,1) double {mustBeInteger, mustBePositive} = 4;
    options.SHT_maxsample (1,1) double {mustBeInteger, mustBePositive} = inf;

    options.SHT_center (:,3) double = [];
    options.SHT_wavnum (1,:) double {mustBeNonnegative} = [];

end

[N_C, N, M] = size(C);
L_C = sqrt(N_C) - 1;
assert(L_C == floor(L_C), 'Invalid size C');

N_D = size(D_pdf, 1);
L_D = sqrt(N_D) - 1;
assert(L_D == floor(L_D), 'Invalid size D_pdf');
assert(size(D_pdf, 2) == 1 || size(D_pdf, 2) == M, 'size(D_pdf, 2) must be 1 or M');

assert(sh_pdf_check(D_pdf, is_real), 'D_pdf is invalid density');

K = complex(zeros([N, N, M]));

if strcmp(options.mode, 'CGC')
    
    D_pdf_res = sh_resize(D_pdf, 2 * L_C); % [(2 * L_C + 1)^2 x 1] or  % [(2 * L_C + 1)^2 x M]

    for m = 1:M     % Iterate over functions
        for i = 1:N % Iterate over channels

            [~, A_i] = sh_mul(sh_conj(C(:, i, m)), C(:, i, m) ); %A: [(2 * L_C + 1)^2 x (L_C + 1)^2]

            E_i_m = A_i * C(:, i:N , m); % [(2 * L_C + 1)^2 x (N - i + 1)]

            if size(D_pdf_res, 2) == 1 % Constant density over M functions
                K(i, i:N, m) = conj(D_pdf_res' * E_i_m); % Conjugate to account for input ordering
            else  % Varied density over M functions
                K(i, i:N, m) = conj(D_pdf_res(:, m)' * E_i_m); % Conjugate to account for input ordering
            end
            K(i:N, i, m) = conj(K(i, i:N, m));

        end

        if options.center_mean
            if size(D_pdf_res, 2) == 1 % Constant density over M functions
                mu = D_pdf(1:min(N_C, N_D))' * C(1:min(N_C, N_D), :, m); % [1 x N]
            else  % Varied density over M functions
                mu = D_pdf(1:min(N_C, N_D), m)' * C(1:min(N_C, N_D), :, m); % [1 x N]
            end
            K(:, :, m) = K(:, :, m) - conj(mu' * mu);
        end
    end
        
elseif strcmp(options.mode, 'SHT')

    assert( ~xor(isempty(options.SHT_center), isempty(options.SHT_wavnum) ), 'options.SHT_center and options.SHT_wavnum must be both empty or both non-empty');

    % Compute number of sampling directions
    if isempty(options.SHT_wavnum) % Minimum twice max_odr due to SH product

        N_C2 = min( options.SHT_maxsample, options.SHT_oversample * (2 * L_C + 1)^2 );

    else % Translation requires higher max_odr expansion

        L_T = floor(max(options.SHT_wavnum) * max(sqrt(sum(options.SHT_center.^2, 2)))  ); % Theoretical maximum expansion order
        N_C2 = min( options.SHT_maxsample, options.SHT_oversample * (2 * L_C + L_T + 1)^2 );

    end
    [theta, phi] = sh_grd_fib(N_C2);

    f_pdf = sh_dec(D_pdf, theta, phi, is_real) / N_C2 * (4 * pi); % [N_C2 x 1] or %[N_C2 x M]

    % Decode expansions along uniform grid
    if isempty(options.SHT_center)

        f_C_full = reshape(sh_dec(reshape(C, [N_C, N * M]), theta, phi, is_real), [N_C2, N, M]); % [N_C2 x N x M] 

    else

        center_wavnum = zeros([N, M, 4]);
        for n = 1:N
            for m = 1:M
                center_wavnum(n, m, :) = [options.SHT_center(n, :), options.SHT_wavnum(m)];
            end
        end
        f_C_full = reshape(sh_dec(reshape(C, [N_C, N * M]), theta, phi, is_real, 'eval', reshape(center_wavnum, [N * M, 4] )), [N_C2, N, M]); % [N_C2 x N x M]            

    end

    for m = 1:M % Iterate over functions

        f_C = f_C_full(:, :, m);

        if size(f_pdf, 2) == 1 % Constant density over M functions
            f_C_pdf = bsxfun(@times, f_C, f_pdf); % [N_C2 x N] 
        else % Varied density over M functions
            f_C_pdf = bsxfun(@times, f_C, f_pdf(:, m)); % [N_C2 x N] 
        end
        K(:, :, m) = conj(f_C' * f_C_pdf); % Conjugate to account for input ordering

        if options.center_mean            
            mu = sum(f_C_pdf, 1); % [1 x N]
            K(:, :, m) = K(:, :, m) - conj(mu' * mu);
        end
    end

else
    error('Unknown options.mode');
end
