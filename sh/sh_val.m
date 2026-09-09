function Y = sh_val(max_odr, theta, phi, is_real, mode)
%Evaluate spherical harmonic (SH) basis functions at spherical coordinates

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%max_odr:       Max SH order 

%theta:         [N x 1]  Co-latitude [0, pi]
%phi:           [N x 1]  Azimuth [0, 2 * pi)

%is_real:       Logical, if true, evaluate real SH

%mode:          String, compute method {'recur', 'direct'}
%                       'recur':    Associated Legendre function recurrence relations
%                       'direct':   Direct evaluation for low-orders (L <= 3),
%                                   otherwise switch to 'recur'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%Y:         [N x (max_odr + 1)^2]   SH evaluated at each spherical coordinate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_val', '-o', 'sh/sh_val_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare recurrence versus direct evaluations

% max_odr = 3;
% rng(6451);
% N = 200;
% theta = rand(N, 1) * pi;
% phi   = rand(N, 1) * 2 * pi;

% Y_recur_cpx       = sh_val(max_odr, theta, phi, false, 'recur');
% Y_direct_cpx      = sh_val(max_odr, theta, phi, false, 'direct');
% err_cpx           = norm(Y_recur_cpx - Y_direct_cpx)

% Y_recur_real      = sh_val(max_odr, theta, phi,  true, 'recur');
% Y_direct_real     = sh_val(max_odr, theta, phi,  true, 'direct');
% err_real          = norm(Y_recur_real - Y_direct_real)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare Matlab versus Mex

% max_odr = 7;
% rng(6451);
% N = 200;
% theta = rand(N, 1) * pi;
% phi   = rand(N, 1) * 2 * pi;

%tic; Y_recur_cpx_mat   = sh_val(max_odr, theta, phi, false, 'recur'); toc_matlab = toc
%tic; Y_recur_cpx_mex   = sh_val_mex(max_odr, theta, phi, false, 'recur'); toc_mex = toc
%err_cpx = norm(Y_recur_cpx_mat - Y_recur_cpx_mex)

%tic; Y_recur_real_mat   = sh_val(max_odr, theta, phi, true, 'recur'); toc_matlab = toc
%tic; Y_recur_real_mex   = sh_val_mex(max_odr, theta, phi, true, 'recur'); toc_mex = toc
%err_real = norm(Y_recur_real_mat - Y_recur_real_mex)

arguments
    max_odr (1,1) double {mustBeNonnegative, mustBeInteger} = 1;

    theta (:,1) double = [0];
    phi   (:,1) double = [0];

    is_real (1,1) logical = false;
    mode (1,:) char {mustBeMember(mode, {'recur', 'direct'})} = 'recur';
end

assert(numel(theta) == numel(phi), 'size theta, phi mismatch');
N = numel(theta);

%Force column vector
theta   = theta(:);
phi     = phi(:);

%Force domain bounds on co-latitude
idx_theta_out = mod(theta, 2 * pi) > pi;
phi(idx_theta_out) = phi(idx_theta_out) + pi;

%Force domain bounds on azimuth
phi = mod(phi, 2 * pi);

%Initialize
Y = complex(zeros(N, (max_odr + 1)^2) );

%Evaluate
for L = 0:max_odr %Iterate over basis orders
    
    if strcmp(mode, 'direct') && L <= 3 %Analytic solutions for low orders
        
        for M = -L:L %Iterate over degrees
            idx = L^2 + M + L + 1;
            Y(:, idx) = sh_val_direct(L, M, theta, phi, is_real);
        end

    else %Use recurrence relation for associated legendre functions

        %Fully normalized associated Legendre removes Condon-Shortley phase  M > 0
        P = asc_legendre(L, cos(theta), 'norm');    
        
        for M = -L:L %Iterate over degrees per order
            
            M_abs = abs(M);
    
            %Compute normalization term
            c = sqrt(1 / (2 * pi) );

            if is_real  %Real
                if M ~= 0
                    c = c * sqrt(2);    %Omit Condon-Shortley phase
                end       

                if M < 0
                    az_comp = sin(M_abs * phi);
                else
                    az_comp = cos(M_abs * phi);                
                end

            else                %Complex
                if M > 0
                    c = c * (-1)^M;     %Add back Condon-Shortley phase for M > 0
                end

                az_comp = exp(1i * M * phi);
            end
        
            idx = L^2 + M + L + 1; %Compute linear index in Y
            Y(:, idx) = c *  (P(M_abs + 1, :).') .* az_comp; %Sample P only at M > 0
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Direct spherical harmonic evaluation upto 3rd order at spherical coordinates

%Input 
%L:             Order
%M:             Degree,  -L <= M <= L
%theta:         [N x 1]  co-latitude (0, pi)
%phi:           [N x 1]  azimuth (0, 2 * pi)

%is_real:       Logical, if true, evaluate real SH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%Y:                     [N x 1]   SH evaluated at each spherical coordinate

function Y = sh_val_direct(L, M, theta, phi, is_real)

arguments
    L = 0;
    M = 0;

    theta = [0];
    phi = [0];

    is_real = false;
end

N = numel(theta);

%Force column vector
theta   = theta(:);
phi     = phi(:);

%Force domain bounds on co-latitude
idx_theta_out = mod(theta, 2 * pi) > pi;
phi(idx_theta_out) = phi(idx_theta_out) + pi;

%Force domain bounds on azimuth
phi = mod(phi, 2 * pi);

%Evaluate
if ~is_real       %Complex
    if L == 0
        Y = 1/2 * sqrt(1 / pi) * ones(N, 1);
    elseif L == 1
        if M == -1
            Y = 1/2 * sqrt(3 / (2 * pi)) * sin(theta) .* exp(-1i * phi);
        elseif M == 0
            Y = 1/2 * sqrt(3 / pi) .* cos(theta);
        elseif M == 1
            Y = -1/2 * sqrt(3 / (2 * pi)) * sin(theta) .* exp(1i * phi);
        else
            Y = nan(N, 1);
            error('Unsupported M');
        end
    elseif L == 2
        if M == -2
            Y = 1/4 * sqrt(15 / (2 * pi)) * sin(theta).^2 .* exp(-2 * 1i * phi);
        elseif  M == -1
            Y = 1/2 * sqrt(15 / (2 * pi)) * sin(theta) .* cos(theta) .* exp(-1i * phi);
        elseif  M == 0
            Y = 1/4 * sqrt(5 / pi) * (3 * cos(theta).^2 - 1);
        elseif  M == 1
            Y = -1/2 * sqrt(15 / (2 * pi)) * sin(theta) .* cos(theta) .* exp(1i * phi);
        elseif  M == 2
            Y = 1/4 * sqrt(15 / (2 * pi)) * sin(theta).^2 .* exp(2 * 1i * phi);
        else
            Y = nan(N, 1);
            error('Unsupported M');
        end
    elseif L == 3
        if M == -3
            Y = 1/8 * sqrt(35 / pi) * sin(theta).^3 .* exp(-3 * 1i * phi);
        elseif  M == -2
            Y = 1/4 * sqrt(105 / (2 * pi)) * sin(theta).^2 .* cos(theta) .* exp(-2 * 1i * phi);
        elseif  M == -1
            Y = 1/8 * sqrt(21 / pi) * sin(theta) .* (5 * cos(theta).^2 - 1) .* exp(-1i * phi);
        elseif  M == 0
            Y = 1/4 * sqrt(7 / pi) * (5 * cos(theta).^3 - 3 * cos(theta));
        elseif  M == 1
            Y = -1/8 * sqrt(21 / pi) * sin(theta) .* (5 * cos(theta).^2 - 1) .* exp(1i * phi);
        elseif  M == 2
            Y = 1/4 * sqrt(105 / (2 * pi)) * sin(theta).^2 .* cos(theta) .* exp(2 * 1i * phi);
        elseif  M == 3
            Y = -1/8 * sqrt(35 / pi) * sin(theta).^3 .* exp(3 * 1i * phi);
        else
            Y = nan(N, 1);
            error('Unsupported M');
        end
    else
        Y = nan(N, 1);
        error('Unsupported L');
    end

else                %Real

    x = sin(theta) .* cos(phi);
    y = sin(theta) .* sin(phi);
    z = cos(theta);
    
    if L == 0
        Y = 1/2 * sqrt(1 / pi);
    elseif L == 1
        if M == -1
            Y = sqrt(3 / (4 * pi)) * y;
        elseif M == 0
            Y = sqrt(3 / (4 * pi)) * z;
        elseif M == 1
            Y = sqrt(3 / (4 * pi)) * x;
        else
            Y = nan(N, 1);
            error('Unsupported M');
        end
    elseif L == 2
        if M == -2
            Y = 1/2 * sqrt(15 / pi) * x .* y;
        elseif M == -1
            Y = 1/2 * sqrt(15 / pi) * y .* z;
        elseif M == 0
            Y = 1/4 * sqrt(5  / pi) * (3 * z.^2 - 1);
        elseif M == 1
        	Y = 1/2 * sqrt(15 / pi) * z .* x;
        elseif M == 2
            Y = 1/4 * sqrt(15 / pi) * (x.^2 - y.^2);
        else
            Y = nan(N, 1);
            error('Unsupported M');
        end
    elseif L == 3
         if M == -3
             Y = 1/4 * sqrt(35 / (2 * pi)) * (3 * x.^2 - y.^2) .* y;
         elseif M == -2
             Y = 1/2 * sqrt(105 / pi) * x .* y .* z;
         elseif M == -1
             Y = 1/4 * sqrt(21 / (2 * pi)) * (5 * z.^2 - 1) .* y;
         elseif M == 0
             Y = 1/4 * sqrt(7 / pi) * (5*z.^3 - 3*z);
         elseif M == 1
             Y = 1/4 * sqrt(21 / (2 * pi)) * (5 * z.^2 - 1) .* x;
         elseif M == 2
             Y = 1/4 * sqrt(105 / pi) * (x.^2 - y.^2) .* z;
         elseif M == 3
             Y = 1/4 * sqrt(35 / (2 * pi)) * (x.^2 - 3 * y.^2) .* x;
         else
            Y = nan(N, 1);
             error('Unsupported M');
         end
    else
        Y = nan(N, 1);
        error('Unsupported L');
    end
end





