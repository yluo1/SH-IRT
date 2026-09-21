function [theta, phi] = sh_pdf_sample(C_pdf, num_samples, mode, is_real)
%Sample spherical coordinates on unit sphere from spherical harmonic expansion
%of probability density function

%Associated Legendre polynomial integrals
%https://functions.wolfram.com/Polynomials/LegendreP2/21/ShowAll.html

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%C:                 [(P + 1)^2 x M] SH coefficients (max order P of M number of functions)
%num_samples:       Number of samples per function, positive
%mode:              String, sampling method {'acceptRejectUniform', 'inverseTransform'}
%                           'acceptRejectUniform':  Accept / Reject from sampling uniform density 
%                                                   over SH weighted upper bounds
%                           'inverseTransform':     Sample theta from inverse transform CDF of co-latitude marginalized over azimuth
%                                                   Sample phi from CDF of azimuth given sampled co-latitude
%is_real:           Logical, if true, evaluate real SH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%theta:             [num_samples x M] Co-latitude (radians)
%phi:               [num_samples x M] Azimuth (radians)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('sh_pdf_sample', '-o', 'sh/sh_pdf_sample_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate bimodal SqExp distribution over sphere and sample spherical coordinates

% rseed = 234;
% max_odr = 5;
% theta = pi/2;
% phi = 0;
% ell = 0.5;
% is_real = false;
% 
% C_pdf = sh_nrm(sh_msq(sh_enc_rbf('SqExp', max_odr, theta, phi, ell, is_real) ...
%                     + sh_enc_rbf('SqExp', max_odr, theta, phi + pi, ell, is_real), is_real), 'Sum');
% sh_plt(C_pdf, 'mercator', is_real);
% 
% num_samples = 1000;

% rng(rseed);
% tic;
% [theta_acceptRejectUniform, phi_acceptRejectUniform] = sh_pdf_sample(C_pdf, num_samples, 'acceptRejectUniform', is_real);
% duration_acceptRejectUniform = toc

% rng(rseed);
% tic;
% [theta_inverseTransform, phi_inverseTransform] = sh_pdf_sample(C_pdf, num_samples, 'inverseTransform', is_real);
% duration_inverseTransform = toc

% sc_plt(theta_acceptRejectUniform(:, 1), phi_acceptRejectUniform(:, 1), ones(num_samples, 1));
% sc_plt(theta_inverseTransform(:, 1), phi_inverseTransform(:, 1), ones(num_samples, 1));

% rng(rseed);
% tic;
% [theta_acceptRejectUniform_mex, phi_acceptRejectUniform_mex] = sh_pdf_sample_mex(C_pdf, num_samples, 'acceptRejectUniform', is_real);
% duration_acceptRejectUniform_mex = toc

% rng(rseed);
% tic;
% [theta_inverseTransform_mex, phi_inverseTransform_mex] = sh_pdf_sample_mex(C_pdf, num_samples, 'inverseTransform', is_real);
% duration_inverseTransform_mex = toc

% err_acceptRejectUniform = norm([theta_acceptRejectUniform - theta_acceptRejectUniform_mex, ...
%                                    phi_acceptRejectUniform - phi_acceptRejectUniform_mex] )

% err_inverseTransform = norm([theta_inverseTransform - theta_inverseTransform_mex, ...
%                              phi_inverseTransform - phi_inverseTransform_mex] )


arguments
     C_pdf (:,:) double {coder.mustBeComplex} = complex(0);
     num_samples (1,1) double {mustBePositive} = 1;
     mode (1,:) char {mustBeMember(mode, {'acceptRejectUniform', 'inverseTransform'})} = 'acceptRejectUniform';

     is_real (1,1) logical = false;
end

[P, M] = size(C_pdf);
P = sqrt(P) - 1;
assert(P == floor(P), 'Invalid size C');

assert(sh_pdf_check(C_pdf, is_real), 'C_pdf expansion not a density');

theta = zeros([num_samples, M]);
phi = zeros([num_samples, M]);

if strcmp(mode, 'acceptRejectUniform')
        
    %Upper bound the PDF via sum of weighted upper bounds of each spherical harmonic basis function
    %sqrt((l-m)! / (l+m)!) *  |P_l^m(x)| <= |P_l(x)|, where  max(|P_l(x)|) = 1 for the associated Legendre Polynomials
    %=> |Y_l^m(theta, phi)| <= sqrt((2*l + 1)/(4*pi)) 
    
    l_vec = zeros(1, (P+1)^2);
    idx = 1;
    for l = 0:P
        for m = -l:l
            l_vec(idx) = l;
            idx = idx + 1;
        end
    end
    Y_ub = sqrt((2 * l_vec + 1) / (4*pi)); %[1 x (P+1)^2] upper bound per SH basis
    ub = Y_ub * abs(C_pdf); %Weight upper bound by modulus of expansion coefficients

    total_samples_drawn = 0;
    for m = 1:M     %Iterate over functions
        for n = 1:num_samples %Iterate over samples            
            while 1

                total_samples_drawn = total_samples_drawn + 1;

                %Sample theta, phi from uniform distribution over sphere
                theta_mn = acos(2 * rand(1) - 1);
                phi_mn = rand(1) * 2 * pi;        
        
                %Evaluate SH
                Y_mn = sh_dec(C_pdf, theta_mn, phi_mn, is_real);
                % if abs(Y_mn) > ub(m)
                %     disp(['invalid: ', num2str(abs(Y_mn) ), ' > ', num2str(ub(m)) ]); %Incorrect upper bound
                % else
                %     %disp(['valid: ', num2str(abs(Y_mn) ), ' <= ', num2str(ub(m)) ]);
                % end

                %Sample from uniform distribution and accept/reject
                if rand(1) <= abs(Y_mn) / ub(m) %Accept
                    theta(n, m) = theta_mn;
                    phi(n, m) = phi_mn;
                    break;
                end
            end
        end
    end

    %disp(['Sampling efficiency: ', num2str((M * num_samples) / total_samples_drawn * 100), ' Percent'] );

elseif strcmp(mode, 'inverseTransform')

    if is_real
        %Convert to complex for CDF sample functions
        C_pdf = sh_re2cpx(C_pdf);
    end

    %Sample theta
    u = rand(num_samples, 1);
    theta = sh_cdf_inv_theta(C_pdf, u, 'mode', 'NewtonRaphsonBisection'); %[num_samples x M]

    %Sample phi
    u = rand(num_samples, 1);
    phi = zeros(num_samples, M);
    for m = 1:M
        phi(:, m) = sh_cdf_inv_phi_cond(C_pdf, u, theta(:, m) , 'mode', 'NewtonRaphsonBisection');
    end
   
else
    error('supported mode');
end

