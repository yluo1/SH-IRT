function [E_pdf] = sh_pdf_transport(C_pdf, D_pdf, t, mode, is_real, options)
%Compute optimal transport displacement interpolations at normalized time t
%Wasserstein metrics between C_pdf to D_pdf

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%C_pdf:     [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%D_pdf:     [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)

%t:         [1 x T] Transport normalized time between [0, 1]

%mode:      String, compute method {'EarthMoverDist', 'SlicedWass', 'LinearInterp', 'GeometricInterp'}
%               'EarthMoverDist':       Minimize earth-mover distance with local non-negative SH expansions over discrete points on sphere.
%                                       Number of local expansions is options.EMD_N_val.
%                                       Solutions are the Monge-Kantorovich displacements of non-negative expansions.
%               'SlicedWass':           Minimize least-squares reconstruction error of SH coefficients under 
%                                       sliced Wasserstein projections of CDFs onto uniform distributed great circles
%               'LinearInterp':         Linear interpolation:       (1-t) * (Y(theta, phi) * C_pdf) + t * (Y(theta, phi) * D_pdf)
%               'GeometricInterp':      Geometric interpolation:    (Y(theta, phi) * C_pdf)^(1-t) * (Y(theta, phi) * D_pdf)^t

%is_real:   Logical, if true, evaluate real SH

%options.EMD_N_val:     Number of discretization points. If 0, then defaults to (P + 1)^2
%options.EMD_K_pdf:     [(P + 1)^2 x 1] Local SH expansion centered on theta = 0, 
%                       [] defaults to squared magnitude projection

%options.SW_fit_mode:   String, fitting method of SH coefficients after sliced Wasserstein projections {'LeastSquares', 'NNLS_MSQ'}
%                           'LeastSquares': Least-squares estimate of SH coefficients
%                           'LS_MSQ':       Least-squares estimate of SH coefficients and fit magnitude-square SH expansion
%                           'NNLS':         Estimate non-negative points on sphere and fit SH expansion
%                           'NNLS_MSQ':     Estimate non-negative points on sphere and fit magnitude-square SH expansion
%options.SW_msq_mode:   String, magnitude squared fitting method for NNLS estimated points of sliced Wasserstein projections {'fmincon', 'sdp'}
%                           See sh_fit_msq.m
%options.SW_N_C_fac:    Oversample number of projections by this non-negative factor
%options.SW_N_u:        Number of samples along marginal CDF

%options.SW_tol:        Inverse CDF tolerance for convergence
%options.SW_max_iter:   Inverse CDF maximum iterations

%options.SW_wt_mode:    String, weighting method per slice {'uniform', 'exp'}
%                           'uniform':      1
%                           'exp':          exp(sh_cdf_theta_diff_sq(R(\Omega) * C_pdf, R(\Omega) * D_pdf))
%                                           Nguyen, K., & Ho, N. (2023). Energy-based sliced wasserstein distance. Advances in Neural Information Processing Systems, 36, 18046-18075.


%options.GI_msq_mode:   String, magnitude squared fitting method {'fmincon', 'sdp'}
%                           See sh_fit_msq.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%E_pdf:     [(P + 1)^2 x M x T] SH coefficients 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Transport sample density functions over sphere

% max_odr = 4;
% theta = pi/2;
% phi = 0;
% ell = 0.5;
% is_real = false;
% t = 0.5;
% EMD_N_val = 200; %Number of discretization points
% %EMD_N_val = 0;
% dB_lim = [-60, 0];
% fig_size = [600, 300];
% disp_phase = false;

% rng(28341);

% % C_pdf = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, theta, phi, ell, is_real)), 'Sum');
% % D_pdf = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, theta, phi + pi/2, ell, is_real)), 'Sum');

% C_pdf = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, theta-pi/4, phi, ell, is_real)), 'Sum');
% D_pdf = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, theta+pi/4, phi, ell, is_real)), 'Sum');

% %C_pdf  = sh_nrm(sh_msq(sh_rand(max_odr, 1, is_real, true)), 'Sum');   %Random function
% %D_pdf  = sh_nrm(sh_msq(sh_rand(max_odr, 1, is_real, true)), 'Sum');   %Random function

% %E_pdf = sh_pdf_transport(C_pdf, D_pdf, t, 'EarthMoverDist', is_real, 'EMD_N_val', EMD_N_val);
% E_pdf = sh_pdf_transport(C_pdf, D_pdf, t, 'SlicedWass', is_real);

% sh_plt(C_pdf, 'mercator', is_real, 'dB_lim', dB_lim, 'disp_phase', disp_phase, 'fig_size', fig_size, 'title_name', 'PDF t = 0');
% sh_plt(E_pdf, 'mercator', is_real, 'dB_lim', dB_lim, 'disp_phase', disp_phase, 'fig_size', fig_size, 'title_name', ['PDF t = ', num2str(t)]);
% sh_plt(D_pdf, 'mercator', is_real, 'dB_lim', dB_lim, 'disp_phase', disp_phase, 'fig_size', fig_size, 'title_name', 'PDF t = 1');

arguments
     C_pdf (:,:) double {coder.mustBeComplex} = complex(0);
     D_pdf (:,:) double {coder.mustBeComplex} = complex(0);
     t (1,:) double {mustBeNonnegative} = 0.5;
     mode (1,:) char {mustBeMember(mode, {'EarthMoverDist', 'SlicedWass', 'LinearInterp', 'GeometricInterp'})} = 'SlicedWass';
     is_real (1,1) logical = false;
     
     options.EMD_N_val (1,1) double {mustBeNonnegative} = 0;
     options.EMD_K_pdf (:,1) double {coder.mustBeComplex} = [];

     options.SW_fit_mode (1,:) char {mustBeMember(options.SW_fit_mode, {'LeastSquares', 'LS_MSq', 'NNLS', 'NNLS_MSq'})} = 'NNLS_MSq';
     options.SW_msq_mode (1,:) char {mustBeMember(options.SW_msq_mode, {'fmincon', 'sdp'})} = 'sdp';
     options.SW_N_C_fac (1,1) double {mustBeNonnegative, mustBeInteger} = 8;
     options.SW_N_u (1,1) double {mustBeNonnegative, mustBeInteger} = 8;
     options.SW_tol (1,1) double {mustBeNonnegative} = 1e-8;
     options.SW_max_iter (1,1) double {mustBeNonnegative} = 1000;
     options.SW_wt_mode (1,:) char {mustBeMember(options.SW_wt_mode, {'uniform', 'exp'})} = 'uniform';

     options.GI_msq_mode (1,:) char {mustBeMember(options.GI_msq_mode, {'fmincon', 'sdp'})} = 'sdp';
     
end

assert(sh_pdf_check(C_pdf, is_real), 'C_pdf is invalid density');
assert(sh_pdf_check(D_pdf, is_real), 'D_pdf is invalid density');
assert(isequal(size(C_pdf), size(D_pdf)), 'C_pdf and D_pdf size mismatch');
assert(all(t <= 1), 't must be between 0 and 1')

[N_C, M] = size(C_pdf);
P = sqrt(N_C) - 1;
assert(P - floor(P) == 0, 'C_pdf invalid size');

T = numel(t);

%Output
E_pdf = zeros([N_C, M, T]);

if strcmp(mode, 'EarthMoverDist') %Earth-mover distance

    %Number of discretization  points
    if options.EMD_N_val == 0
        N = N_C; 
    else
        N = options.EMD_N_val;
    end

    if ~isempty(options.EMD_K_pdf)
        assert(size(options.EMD_K_pdf, 1) == N_C, 'Size EMD_K_pdf mismatch C_pdf');
    end

    %Generate uniform points over sphere
    [theta, phi] = sh_fib(N);
    v = zeros([3, N]);
    [v(1,:), v(2,:), v(3,:)] = sph2cart(phi, pi/2 - theta, ones(N, 1));
    
    %Compute geodesic (angle) distance between discretized spherical coordinates
    D = acos(max(min(v'*v, 1), -1)); %[N x N] [0 to pi]

    %Decopose C_pdf, D_pdf into sum of local pdf kernels
    if isempty(options.EMD_K_pdf)
        C_kers_real = sh_enc_proj_msq(P, theta, phi, true);     %[(max_odr + 1)^2 x N], local kernel expansions at each spherical coordinates, real        
    else
        C_kers_real = zeros([N_C, N]);
        for n = 1:N
            C_kers_real(:, n) = sh_rot(options.EMD_K_pdf, [theta(n), phi(n), 0], true);
        end
    end
    %Integrated least-squares, convert SH to real form
    w_C = zeros([N, M]);
    w_D = zeros([N, M]); 
    if ~is_real
        C_pdf_real = real(sh_cpx2re(C_pdf));
        D_pdf_real = real(sh_cpx2re(D_pdf));
    else
        C_pdf_real = real(C_pdf);
        D_pdf_real = real(D_pdf);
    end    
    for m = 1:M
        [w_C(:, m), resnorm_C, residual_C] = lsqnonneg(C_kers_real, C_pdf_real(:, m)); %[N x 1] non-negative weights per kernel
        [w_D(:, m), resnorm_D, residual_D] = lsqnonneg(C_kers_real, D_pdf_real(:, m)); %[N x 1] non-negative weights per kernel
    end

    % sh_plt(C_pdf, 'mercator', is_real, 'dB_lim', [-60, 0]);
    % sh_plt(C_kers_real * w_C, 'mercator', true, 'dB_lim', [-60, 0]);
    ;

    %Assignment: X = [N x N], moves X(i,j) mass from 
    %C_pdf expansion at [theta(i), phi(i)] to 
    %D_pdf expansion at [theta(j), phi(j)]

    %Iterate over time
    for idx_t = 1:T

        %Determine interpolated unit-vector V_t(i,j)    
        V_t = zeros([3, N, N]);
        for i = 1:N
            for j = 1:N
                if abs(D(i, j)) <= 1e-8
                    V_t(:, i, j) = v(:, i);
                else
                    V_t(:, i, j) = (sin((1 - t(idx_t)) * D(i, j)) * v(:, i) + sin(t(idx_t) * D(i, j)) * v(:, j)) / sin(D(i, j));
                    V_t(:, i, j) = V_t(:, i, j) / norm(V_t(:, i, j));  %Normalize to account for any round-off error
                end
            end
        end
        V_t_list = reshape(V_t, [3, N * N]); 
    
        %Convert into spherical coordinates
        [phi_t_list, theta_t_list] = cart2sph(V_t_list(1, :), V_t_list(2, :), V_t_list(3, :)); %[1 x  N * N]
        theta_t_list = pi/2 -  theta_t_list;

        %Iterate over functions
        for m = 1:M

            %Objective
            f = D(:);

            %Inequality constraints
            A_in_D = ones(1, N);
            for n = 2:N
                A_in_D = blkdiag(A_in_D, ones(1,N));
            end
            b_in_D = w_D(:, m);
    
            A_out_C = repmat([1, zeros(1, N-1)], [1, N]);
            for n = 2:N
                A_out_C = [A_out_C; circshift(A_out_C(n-1, :), 1)];
            end
            b_in_C = w_C(:, m);
            
            A = [A_in_D; A_out_C];
            b = [b_in_D; b_in_C];
    
            %Equality constraints for preserving mass, C and D expansion mass converges to 
            %same quantity for pdfs as N -> inf
            Aeq = ones(1, N * N);
            beq = min(sum(w_C(:, m)), sum(w_D(:, m)));
    
            %Bounds
            lb = zeros(N * N, 1);
            ub = beq * ones(N * N, 1);
            
            %Solve
            [x, w_m_idx_t, exitflag] = linprog(f,A,b,Aeq,beq,lb,ub); %[N * N x 1]
            
            %Expand
            if isempty(options.EMD_K_pdf)
                C_kers_t_list = sh_enc_proj_msq(P, theta_t_list, phi_t_list, is_real);
            else
                C_kers_t_list = zeros([N_C, N * N]);
                for n = 1:(N*N)
                    C_kers_t_list(:, n) = sh_rot(options.EMD_K_pdf, [theta_t_list(n), phi_t_list(n), 0], is_real);
                end
            end
            E_pdf(:, m, idx_t) = C_kers_t_list * x;
            E_pdf(:, m, idx_t) = sh_nrm(E_pdf(:, m, idx_t), 'Sum'); %Normalize
    
        end
    end

elseif strcmp(mode, 'SlicedWass') %Sliced Wasserstein projections

    N   = options.SW_N_C_fac * N_C; %Number of projections
    N_u = options.SW_N_u;           %Number of points to sample in CDF

    rot_intrinsic = false;

    [theta, phi] = sh_fib(N);
    %psi = rand(N, 1) * 2 * pi;
    %psi = linspace(0, 2 * pi * (N / N_C), N)';
    psi = zeros(N, 1);

    if ~is_real
        %Convert from complex to real
        C_pdf = sh_cpx2re(C_pdf);
        D_pdf = sh_cpx2re(D_pdf);        
    end

    %CDF targets
    u_list = (1:N_u)' ./ (N_u + 1);
    CDF_u = repmat(u_list, [1, N]);

    for m = 1:M %Iterate over functions
        C_pdf_rot_m = zeros([N_C, N]); %[(P+1)^2 x N]
        D_pdf_rot_m = zeros([N_C, N]);
        A = zeros([N_u, N, N_C, T]);

        wt = ones([N, T]); %Weight per projection 

        %Iterate over projections
        for n = 1:N

            %Compute rotation matrices, real
            C_pdf_rot_m(:, n) = sh_rot(C_pdf(:, m), -[phi(n), theta(n), psi(n)], true, rot_intrinsic);
            D_pdf_rot_m(:, n) = sh_rot(D_pdf(:, m), -[phi(n), theta(n), psi(n)], true, rot_intrinsic);
            [R_list] = sh_rot_mat(P, [phi(n), theta(n), psi(n)], true, 'convention', 'zyx', 'rot_intrinsic', rot_intrinsic);

            %Non-uniform weighting
            if strcmp(options.SW_wt_mode, 'exp')
                W2_n  = sh_cdf_inv_theta_diff_sq(C_pdf_rot_m(:, n), D_pdf_rot_m(:, n));
                wt(n) = exp(W2_n);
            end

            %Iterate over time
            for idx_t = 1:T

                for i = 1:N_u %Iterate over theta(u) = CDF_inv(u, t) at displacement time t

                    %Compute CDF_inv(u, t)
                    theta_CDF = (1-t(idx_t)) * sh_cdf_inv_theta(C_pdf_rot_m(:, n), u_list(i), 'tol', options.SW_tol, 'max_iter', options.SW_max_iter) ...
                              + t(idx_t)     * sh_cdf_inv_theta(D_pdf_rot_m(:, n), u_list(i), 'tol', options.SW_tol, 'max_iter', options.SW_max_iter);
                    cos_theta = cos(theta_CDF);
                  
                    for l = 0:P %Rows are coefficients for computing CDF
                        idx = (l + 1)^2 - l;

                        if l == 0
                            A(i, n, idx, idx_t) = sqrt(pi / (2*l + 1) ) * ( 1 - legendreP2(l+1, cos_theta)  );
                        else
                            A(i, n, idx, idx_t) = sqrt(pi / (2*l + 1) ) * ( (legendreP2(l-1, cos_theta) - legendreP2(l+1, cos_theta) ) );
                        end

                        idx_rot = ((l^2) + 1):((l+1)^2);
                        A(i, n, idx_rot, idx_t) = squeeze(A(i, n, idx_rot, idx_t)).' * R_list(l+1).R;                        
                    end
                    
                end
                
            end

        end
            
        %Iterate over time
        for idx_t = 1:T
            A_mat = reshape(A(:, :, :, idx_t), [N_u * N, N_C]);
            cond_A = cond(A_mat)

            %Setup weighted least squares
            wt_idx_t = reshape(ones(N_u, 1) * wt(:, idx_t)', [N_u * N, 1]); %[N_u * N x 1]
            wt_idx_t = wt_idx_t / sum(wt_idx_t); %Normalize weights
            wA_mat = bsxfun(@times, A_mat, wt_idx_t);
            wCDF_u = CDF_u(:) .* wt_idx_t;
            
            if strcmp(options.SW_fit_mode, 'LeastSquares')   %Least squares fit, does not ensure non-negative pdf
                E_pdf(:, m, idx_t) = wA_mat \ wCDF_u(:);        
                   
            elseif strcmp(options.SW_fit_mode, 'LS_MSq') %Least-squares fit, and magnitude squared fit
              
                [theta_fit, phi_fit] = sh_fib(N_C);
                E_pdf_LS = wA_mat \ wCDF_u(:); 
                E_pdf_LS = sh_nrm(E_pdf_LS, 'Sum'); %Normalize, ~ 1/(2*pi)
                X = real(sh_dec(E_pdf_LS, theta_fit, phi_fit, true));
                if strcmp(options.SW_msq_mode, 'fmincon')
                    C0 = sh_rand(floor(P/2), 100, true);
                else
                    C0 = [];
                end
                [E_pdf(:, m, idx_t), ~, err] = sh_fit_msq(X, theta_fit, phi_fit, P, true, options.SW_msq_mode, 'C0', C0); 
                
                err

            elseif strcmp(options.SW_fit_mode, 'NNLS') %NNLS estimated points on sphere
              
                [theta_fit, phi_fit] = sh_fib(N_C);
                Y_fit = sh_val(P, theta_fit, phi_fit, true);  %[N_C x N_C]        
                A_inv_Y_fit = wA_mat / Y_fit;
                X = lsqnonneg(A_inv_Y_fit, wCDF_u(:));    %min ||A_mat * inv(Y) * X -  CDF_u||^2, s.t. x >= 0
                E_pdf(:, m, idx_t) = Y_fit \ X;
    
            elseif strcmp(options.SW_fit_mode, 'NNLS_MSq') %NNLS estimated points on sphere, and magnitude squared fit
              
                [theta_fit, phi_fit] = sh_fib(N_C);
                Y_fit = sh_val(P, theta_fit, phi_fit, true);  %[N_C x N_C]        
                A_inv_Y_fit = wA_mat / Y_fit;
                X = lsqnonneg(A_inv_Y_fit, wCDF_u(:));    %min ||A_mat * inv(Y) * X -  CDF_u||^2, s.t. x >= 0
                %X = X / (2*pi); %Normalize
                if strcmp(options.SW_msq_mode, 'fmincon')
                    C0 = sh_rand(floor(P/2), 100, true);
                else
                    C0 = [];
                end
                [E_pdf(:, m, idx_t), ~, err] = sh_fit_msq(X, theta_fit, phi_fit, P, true, options.SW_msq_mode, 'C0', C0); 
                
                err

            else
                error('Unsupported options.SW_mode');
            end
            E_pdf(:, m, idx_t) = sh_nrm(E_pdf(:, m, idx_t), 'Sum'); %Normalize
            
            if ~is_real %Convert to complex expansions
                E_pdf(:, m, idx_t) = sh_re2cpx(E_pdf(:, m, idx_t));
            end
        end
    end

elseif strcmp(mode, 'LinearInterp')

    %Iterate over time
    for idx_t = 1:T
        E_pdf(:, :, idx_t) = (1 - t(idx_t)) * C_pdf + t(idx_t) * D_pdf;
        E_pdf(:, :, idx_t) = sh_nrm(E_pdf(:, :, idx_t), 'Sum');
    end

elseif strcmp(mode, 'GeometricInterp')
    
    [theta, phi] = sh_fib(N_C);
    f_C = real(sh_dec(C_pdf, theta, phi, is_real));
    f_D = real(sh_dec(D_pdf, theta, phi, is_real));
    
    %Iterate over time
    for idx_t = 1:T
        X = (f_C.^(1-t(idx_t))) .* (f_D.^t(idx_t));
        [E_pdf(:, :, idx_t), ~, err] = sh_fit_msq(X, theta, phi, P, is_real, options.GI_msq_mode); 
        E_pdf(:, :, idx_t) = sh_nrm(E_pdf(:, :, idx_t), 'Sum');
    end
 
else
    error('unknown mode');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = legendreP2(n, x)
arguments
    n (1,1) double = 0;
    x (:,1) double = 0;
end
P = asc_legendre(n, x);
Y = P(1,:).';
