function f = srbf_val(mode, theta_s, phi_s, ell, theta, phi)
%Evaluate spherical radial basis functions (RBF) at spherical coordinates, and bandwidth

%Reference:
%Y. Luo, "Spherical harmonic covariance and magnitude function encodings for beamformer design," 
%EURASIP Journal on Audio, Speech, and Music Processing. 2021. 10.1186/s13636-021-00230-7. 

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%mode:          String, RBF function

%               'SqExp':                Squared exponential
%               'SqExpNorm':            Normalized squared exponential
%               'Mat52':                Matern (\nu = 5/2)
%               'Mat32':                Matern (\nu = 3/2)
%               'Exp':                  Exponential
%               'ExpNorm':              Normalized exponential
%               'Sinc':                 Sinc (unnormalized)

%theta_s:       Co-latitude [0, pi] steering direction
%phi_s:         Azimuth center [0, 2 * pi) steering direction
%ell:           Scalar, bandwidth parameter, positive, same units as d

%theta:         [N x 1] Co-latitude [0, pi]
%phi:           [N x 1] Azimuth [0, 2 * pi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%f:             [N x 1] Function evaluations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation

%codegen('srbf_val', '-o', 'sh/srbf_val_mex')


arguments
    mode (1,:) char  = 'SqExp';
    theta_s (1,1) double = 0;
    phi_s (1,1) double = 0;
    ell (1,1) double {mustBePositive} = 1;
    theta (:,1) double  = 0;
    phi (:,1) double  = 0;
end

u = zeros(3, 1);
[u(1), u(2), u(3)] = sph2cart(phi_s, pi/2 - theta_s, 1);

assert(numel(theta) == numel(phi), 'theta, phi size mismatch');

N = numel(theta);
v = zeros(3, N);
[v(1,:), v(2,:), v(3,:)] = sph2cart(phi, pi/2 - theta, ones(size(theta)));

%Compute distance
x = acos(max(min(v' * u, 1), -1) );
d = 2 * sin(abs(x)/2);

if contains(mode, 'SqExp')

    if contains(mode, 'SqExpNorm')
        c = 1 / (2 * pi * ell^2 * (1 - exp(-2 / (ell^2))));
    else
        c = 1;
    end

    f = c * exp(-d.^2 ./ (2 * ell.^2));   

elseif strcmp(mode, 'Mat52')

    f = (1 + sqrt(5) .* d ./ ell + 5/3 .* (d ./ ell).^2 ) .* exp(-sqrt(5) .* d ./ ell);    

elseif strcmp(mode, 'Mat32')

    f = (1 + sqrt(3) .* d ./ ell) .* exp(-sqrt(3) * d / ell );    

elseif contains(mode, 'Exp')

    if contains(mode, 'ExpNorm')
        c =  1 / (2 * pi * (ell^2 - ell * exp(-2 / ell) * (ell + 2)));
    else
        c = 1;
    end

    f = c * exp(-(d ./ ell));
    
elseif contains(mode, 'Sinc')

    f = sinc(ell * d / pi);

else

    f = nan(size(d));
    error('Unknown mode');

end
