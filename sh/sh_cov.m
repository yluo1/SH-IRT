function K = sh_cov(C, D_pdf, is_real, options)
%Compute covariance between spherical harmonic expansions in C over density D_pdf

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:             [(L_C + 1)^2 x N x M]   SH coefficients (max order P of N number of expansions, M functions)
%D_pdf:         [(L_D + 1)^2 x 1        Density SH coefficients
%is_real:       Logical, if true, evaluate real SH

%options.mode:          String, integration (SH multiplication) method {'CGC', 'SHT'}
%                           'CGC':   Clebsch-Gordan coefficients
%                           'SHT':   Spherical harmonic transform (inverse -> prod -> forward)

%options.center_mean:       Logical, if true, subtract mean from covariance matrix
%options.SHT_oversample:    Oversample factor for mode = 'SHT'

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


arguments
    C (:,:) double {coder.mustBeComplex} = complex(0);
    D_pdf (:,:) double {coder.mustBeComplex} = complex(0);
    is_real (1,1) logical = false;
    
    options.mode (1,:) char {mustBeMember(options.mode, {'CGC', 'SHT'})} = 'SHT';
    options.center_mean (1,1) logical = false;

    options.SHT_oversample (1,1) double {mustBeInteger, mustBePositive} = 4;
end

[N_C, N, M] = size(C);
L_C = sqrt(N_C) - 1;
assert(L_C == floor(L_C), 'Invalid size C');

N_D = size(D_pdf, 1);
L_D = sqrt(N_D) - 1;
assert(L_D == floor(L_D), 'Invalid size D');

assert(sh_pdf_check(D_pdf, is_real), 'D_pdf is invalid density');

K = complex(zeros([N, N, M]));

if strcmp(options.mode, 'CGC')

    D_pdf_res = sh_resize(D_pdf, 2 * L_C); % [(2 * L_C + 1)^2 x 1]

    for m = 1:M    
        for i = 1:N

            [~, A_i] = sh_mul(sh_conj(C(:, i, m)), C(:, i, m) ); %A: [(2 * L_C + 1)^2 x (L_C + 1)^2]

            E_i_m = A_i * C(:, i:N , m); % [(2 * L_C + 1)^2 x (N - i + 1)]

            K(i, i:N, m) = conj(D_pdf_res' * E_i_m); % Conjugate to account for input ordering
            K(i:N, i, m) = conj(K(i, i:N, m));

        end

        if options.center_mean
            mu = D_pdf(1:min(N_C, N_D))' * C(1:min(N_C, N_D), :, m); % [1 x N]
            K(:, :, m) = K(:, :, m) - conj(mu' * mu);
        end
    end
        
elseif strcmp(options.mode, 'SHT')
    
    N_C2 = options.SHT_oversample * (2 * L_C + 1)^2;

    [theta, phi] = sh_grd_fib(N_C2);

    f_pdf = sh_dec(D_pdf, theta, phi, is_real) / N_C2 * (4 * pi); %[N_C2 x 1]

    for m = 1:M

        f_C = sh_dec(C(:, :, m), theta, phi, is_real); % [N_C2 x N] 
        f_C_pdf = bsxfun(@times, f_C, f_pdf); % [N_C2 x N] 
        K(:, :, m) = conj(f_C' * f_C_pdf); % Conjugate to account for input ordering
        
        if options.center_mean            
            mu = sum(f_C_pdf, 1); % [1 x N]
            K(:, :, m) = K(:, :, m) - conj(mu' * mu);
        end
    end

else
    error('Unknown options.mode');
end
