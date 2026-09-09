function [C, D, loglike] = sh_pdf_est(theta, phi, max_odr, is_real, mode, options)
%Fit sum of magnitude squared density function spherical harmonic expansion to 
%sampled spherical coordinates (theta, phi)

%f_D(theta, phi) = Y(theta, phi) * D = sum_{n=1}^{N_S} (|Y(theta, phi) * C(:, n)|^2)

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%theta:         [N x 1] Co-latitude [0, pi]
%phi:           [N x 1] Azimuth [0, 2 * pi)

%max_odr:       Max SH order 

%is_real:       Logical, if true, evaluate real SH

%mode:          Fitting method, {'MaxPow', 'MaxPowCosWt', 'MaxPowUniRegu', 'TotalVariance', 'KernelDensityProj', 'MaxLogLike'}
%                               'Uniform'               Uniform distribution
%                               'MaxPow'                Maximum power
%                               'MaxPowCosWt'           Maximum power, cosine direction weighting
%                               'MaxPowUniWt'           Maximum power, weight non-omni-directional component
%                                                       by options.MaxPowUniWt_alpha [0 to 1]
%                               'KernelDensityProjPre'  Kernel density estimate via pre-sum of SH projections into C
%                               'KernelDensityProjPost' Kernel density estimate via post-sum of magnitude squared SH projections into D
%                               'KernelDensitySqProj'   Kernel density estimate via sum of SH magnitude squared projection
%                               'PrincipalAxes'         Eigenvectors weighted by sqrt(eigenvalues)
%                               'SpectralDecomp'        Spectral decomposition, C is eigenvectors weighted by sqrt(eigenvalues)
%                               'MaxLogLike'            Maximum log-likelihood

%options:       struct

%options.MaxLogLike_mode:               String, fitting method {'fmincon', 'sdp'}
%                                               'fmincon':  nonlinear constrained least squares
%                                               'sdp':      Semi-definite program sum of squares (cvx)
%options.MaxLogLike_options_fmincon:    options struct for 'max_loglike' fmincon 
%options.MaxLogLike_x0_list:            [(max_odr + 1)^2 x M] M initial starting points, 
%                                                      otherwise, start with kernel density estimates
%options.MaxLogLike_relax_eq_constr:    Logical, if true, convert unity energy equality constraint to inequality

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%C:         [(max_odr + 1)^2 x N_S] SH coefficients
%D:         [(2 * max_odr + 1)^2 x 1] SH coefficients
%loglike:   Log-likelihood

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate small set of spherical coordinates and fit density functions 

%rng(11213);
%max_odr = 6;
%is_real = false;
%N = 10;
%theta = acos(2 * rand(N, 1) - 1);
%phi = rand(N, 1) * 2 * pi;
%dB_lim = [-60, 0];

% %theta = [pi/2; pi/2; pi/2; pi/2];
% %phi = [0; pi/4; pi/8;  pi];

%theta = acos(2 * rand(N, 1) - 1) / 2; %[0 to pi/2]
%phi = rand(N, 1) * 2 * pi;
%theta = [theta; pi]; %Add south pole point
%phi = [phi; 0];

%mode_list = {'Uniform', 'KernelDensityProjPre', 'KernelDensityProjPost', 'KernelDensitySqProj', 'PrincipalAxes', 'MaxLogLike'};
%N_modes = numel(mode_list);
%C_list = cell(1, N_modes);
%D_list = cell(1, N_modes);
%loglike_list = cell(1, N_modes);
%MaxLogLike_x0_list = randn((max_odr + 1)^2, 50);
%%MaxLogLike_x0_list = [];
%for i = 1:N_modes
    %[C_list{i}, D_list{i}, loglike_list{i}] = sh_pdf_est(theta, phi, max_odr, is_real, mode_list{i}, 'MaxLogLike_x0_list', MaxLogLike_x0_list);
    %sh_plt(D_list{i}, 'mercator', is_real, 'disp_phase', false, 'disp_theta_phi', [theta, phi], 'dB_lim', dB_lim, 'title_name', mode_list{i});
    %loglike_list{i}
%end

% tic;
% [theta_sample_ar, phi_sample_ar] = sh_pdf_sample(D_list{3}, 250, 'acceptRejectUniform', is_real);
% toc
% sh_plt(D_list{3}, 'mercator', is_real, 'disp_theta_phi', [theta_sample_ar, phi_sample_ar], 'disp_theta_phi_markersize', 16, 'disp_phase', false);
% 
% tic;
% [theta_sample_invt, phi_sample_invt] = sh_pdf_sample(D_list{3}, 250, 'inverseTransform', is_real);
% toc
% sh_plt(D_list{3}, 'mercator', is_real, 'disp_theta_phi', [theta_sample_invt, phi_sample_invt], 'disp_theta_phi_markersize', 16, 'disp_phase', false);

arguments
    theta (:,1) double = [0];
    phi   (:,1) double = [0];

    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;

    is_real (1,1) logical = false;
    
    mode (1,:) char {mustBeMember(mode, {'Uniform', 'MaxPow', 'MaxPowCosWt', 'MaxPowUniWt', ...
        'KernelDensityProjPre', 'KernelDensityProjPost', 'KernelDensitySqProj', 'PrincipalAxes', 'SpectralDecomp', 'SpectralMaxLogLike', 'MaxLogLike'})} = 'MaxLogLike';

    options.MaxLogLike_options_fmincon = optimoptions("fmincon", ...
        Algorithm="interior-point", EnableFeasibilityMode=false, ...
        SpecifyObjectiveGradient=true, SpecifyConstraintGradient=true, ...
        Display="off", ScaleProblem=true, ...
        FunctionTolerance=1e-8, ConstraintTolerance=1e-8, OptimalityTolerance=1e-10, StepTolerance=1e-8, ...
        MaxIterations=1000, MaxFunctionEvaluations = 10000);

    options.MaxLogLike_mode (1,:) char {mustBeMember(options.MaxLogLike_mode, {'fmincon', 'sdp'})} = 'fmincon';
    options.MaxLogLike_x0_list (:,:) double = [];
    options.MaxLogLike_relax_eq_constr (1,1) logical = false;

    options.MaxPowUniWt_alpha (1,1) double {mustBeNonnegative} = 0.2424;

end

assert(numel(theta) == numel(phi), 'Size mismatch theta, phi');
N = numel(theta);
N_C = (max_odr + 1)^2;

if strcmp(mode, 'Uniform')
    C = sh_enc_uni(max_odr);
    C = C / sqrt(C'*C);

elseif  strcmp(mode, 'MaxPow')
    Y = sh_val(max_odr, theta, phi, is_real); %[N x (max_odr + 1)^2]
    [C, eig_val] = eigs(Y'*Y, 1, 'largestabs');
   
elseif strcmp(mode, 'MaxPowCosWt')
    Y = sh_val(max_odr, theta, phi, is_real); %[N x (max_odr + 1)^2]
    
    %Cosine similarity weighting between directions
    V = zeros(N, 3);
    [V(:, 1), V(:, 2), V(:, 3)] = sph2cart(phi, pi/2 - theta, ones(N, 1));    
    T = V' * Y; %[3 x (max_odr + 1)^2]    
    [C, eig_val] = eigs(T' * T, 1, 'largestabs');

elseif strcmp(mode, 'MaxPowUniWt')
    assert(options.MaxPowUniWt_alpha <= 1, 'options.MaxPowUniWt_alpha must be between [0, 1]');

    Y = sh_val(max_odr, theta, phi, is_real); %[N x (max_odr + 1)^2]    

    W =  diag([1, options.MaxPowUniWt_alpha * ones(1, N_C - 1) ]);
    Y = Y * W;
    [C, eig_val] = eigs(Y'*Y, 1, 'largestabs');

elseif strcmp(mode, 'KernelDensityProjPre') %C is sum of projections
    C = zeros([N_C, 1]);
    %Initial point: Kernel density estimate from projections onto SH
    for n = 1:N     %Kernel density estimate
        C = C + sh_enc_proj(max_odr, theta(n), phi(n), is_real);
    end
    C = C / sqrt(C'*C);

elseif strcmp(mode, 'KernelDensityProjPost') %C is list of projections
    %Same solution as SpectralDecomp with less computation
    C = sh_enc_proj(max_odr, theta, phi, is_real);  
    C = C / sqrt(sum(sum(conj(C) .* C)));
    ;

elseif strcmp(mode, 'KernelDensitySqProj')
    C = zeros([N_C, 1]);
    %Initial point: Kernel density estimate from magnitude squared projections onto SH
    for n = 1:N     %Kernel density estimate
        C = C + sh_enc_proj_msq(max_odr, theta(n), phi(n), is_real);
    end
    C = C / sqrt(C'*C);

elseif strcmp(mode, 'PrincipalAxes')
     
    Y = sh_val(max_odr, theta, phi, true); %[N x (max_odr + 1)^2], reals
    [C_list, eig_val] = eig(Y'*Y);
    C_list = bsxfun(@times, C_list, sh_int(C_list, 'Sum')); %Resolve sign ambiguity and normalize
    C = C_list * sqrt(abs(diag(eig_val)));
    %sh_plt(C, 'mercator', true)

    if ~is_real
        C = sh_re2cpx(C);
    end
    C = C / sqrt(C'*C);

elseif  strcmp(mode, 'SpectralDecomp')
    Y = sh_val(max_odr, theta, phi, is_real); %[N x (max_odr + 1)^2]
    [eig_vecs, eig_vals] = eig(Y'*Y);
    
    C = eig_vecs * sqrt(eig_vals);
    C = C / sqrt(sum(sum(conj(C) .* C)));

elseif  strcmp(mode, 'SpectralMaxLogLike')
    Y = sh_val(max_odr, theta, phi, is_real); %[N x (max_odr + 1)^2]
    [V, D] = eig(Y'*Y);
    
    cvx_clear
    cvx_begin
        variable w(N_C, 1)
        expression x(N, 1)

        for n = 1:N
            x(n) = trace(w' * diag(V' * (Y(n,:)' * Y(n,:)) * V) );
        end
            
        maximize(sum_log(x))

        subject to
            w'*w <= 1;
            w >= 0;

    cvx_end

    C = V * diag(w);
    ;

elseif strcmp(mode, 'MaxLogLike')

    Y = sh_val(max_odr, theta, phi, true); %[N x (max_odr + 1)^2], Real SH
    
    if strcmp(options.MaxLogLike_mode, 'fmincon')

        if isempty(options.MaxLogLike_x0_list)       
            x0_list = zeros([N_C, 1]);
            %Initial point: Kernel density estimate from projections onto SH
            for n = 1:N     %Kernel density estimate
                x0_list(:, 1) = x0_list(:, 1) + sh_enc_proj(max_odr, theta(n), phi(n), is_real);
            end        
            %Display
            % sh_plt(x0_list(:, 1), 'mercator', is_real, 'disp_theta_phi', [theta, phi], 'disp_phase', false);
     
            if ~is_real
                x0_list = sh_cpx2re(x0_list);
            end
    
        else %Custom starting points
            assert(size(options.MaxLogLike_x0_list, 1) == size(Y, 2), 'optins.MaxLogLike_x0_list size mismatch');
            x0_list = options.MaxLogLike_x0_list;
    
            if ~is_real
                x0_list = real(sh_cpx2re(x0_list));
            end
        end
        %Normalize initial guesses to unit circle
        x0_list = bsxfun(@rdivide, x0_list, sqrt(sum(x0_list .* conj(x0_list), 1)));
       
        %Find maximum likelihood solution from set of initial guesses
        M = size(x0_list, 2);
    
        %Initial guess
        C = x0_list(:, 1);
        max_loglike = -f_neg_log_likelihood(C, sh_val(max_odr, theta, phi, true));
    
    
        for m = 1:M        
            x0 = x0_list(:, m);
            loglike_x0 = -f_neg_log_likelihood(x0, Y);
            if max_loglike < loglike_x0
                max_loglike = loglike_x0;
                C = x0;
            end
          
            [C_m, fval, exitflag]= fmincon(@(x) f_neg_log_likelihood(x, Y), x0, ...
                [],[],[],[],[],[], @(x) sh_msq_pdf_nonlcon(x, options.MaxLogLike_relax_eq_constr), options.MaxLogLike_options_fmincon);
    
            
            C_m = C_m / sqrt(C_m' * C_m);
            loglike_m = -f_neg_log_likelihood(C_m, Y);
    
            if max_loglike < loglike_m
                max_loglike = loglike_m;
                C = C_m;
            end
        end

    elseif strcmp(options.MaxLogLike_mode, 'sdp') %Semi-definite programming relaxation (cvx)
        
        cvx_clear
        cvx_begin
            variable Q(N_C, N_C)
            expression x(N, 1)

            for n = 1:N
                x(n) = trace(Q * (Y(n,:)' * Y(n,:)));
            end
                
            maximize(sum_log(x))

            subject to
                Q == semidefinite(N_C);
                trace(Q) <= 1; 

        cvx_end

        %Check tightness
        [C_eig_vec, C_eig_val] = eigs(double(Q), 1, 'largestabs'); %Largest eigenvector-eigenvalue pair       
        if C_eig_val / trace(Q) <= 0.95
            disp(['Warning: Not tight as largest eigenvalue / trace: ', num2str(C_eig_val / trace(Q))] );
        end

        %Reconstruct D as sum of squares
        [C_eig_vec, C_eig_val] = eig(double(Q)); %Weighted eigenvector-eigenvalue pair
        C_eig_val = diag(C_eig_val);
        idx_eigs = (C_eig_val / trace(Q) > 1e-8);
        C = C_eig_vec(:, idx_eigs) * diag(sqrt(C_eig_val(idx_eigs))); 

    else
        error('Unsupported options.MaxLogLike_mode');
    end

    if ~is_real %Convert back to original bases
        C = sh_re2cpx(C);
    end

else
    error('Unsupported mode');
end

D = sum(sh_msq(C, is_real), 2);

%Compute log-likelihood
if nargout > 2
    Y_loglike = sh_val(max_odr, theta, phi, is_real); %[N x (max_odr + 1)^2]
    loglike = -f_neg_log_likelihood(C, Y_loglike);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Negative log-likelihood and gradient

%Input
%C:             [N_C x N_S]
%Y:             [N x N_C]

%Output
%f:             Scalar, negative log-likelihood
%g:             [N_C x N_S] Gradient

function [f, g] = f_neg_log_likelihood(C, Y)
[N, N_C] = size(Y);
N_S = size(C, 2);
f = 0;
g = zeros(N_C, N_S);
for n = 1:N
    YC_n = Y(n,:) * C; %[1 x N_S]
    %YC_msq_n = trace(YC_n'*YC_n);
    YC_msq_n = sum(conj(YC_n) .* YC_n);
    f = f + log(YC_msq_n);
    g = g + Y(n,:)' * YC_n / YC_msq_n;
end
g = 2 * g;

%Negative log-likelihood
f = -f;
g = -g;
