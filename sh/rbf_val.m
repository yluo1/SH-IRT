function f = rbf_val(mode, d, ell, enable_disp)
%Evaluate radial basis functions (RBF) at distance and bandwidth

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

%d:             [N x M] Distance, non-negative
%ell:           Scalar, bandwidth parameter, positive, same units as d

%enable_disp:   Logical, if true, plot function evaluations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%f:             [N x M] Function evaluations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation

%codegen('rbf_val', '-o', 'sh/rbf_val_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare function types

% d = linspace(0, 10, 256);
% ell = 1;
% f_SqExp        = rbf_val('SqExp', d, ell);
% f_SqExpNorm    = rbf_val('SqExpNorm', d, ell);
% f_Mat52        = rbf_val('Mat52', d, ell);
% f_Mat32        = rbf_val('Mat32', d, ell);
% f_Exp          = rbf_val('Exp', d, ell);
% f_ExpNorm      = rbf_val('ExpNorm', d, ell);
% f_Sinc         = rbf_val('Sinc', d, ell);
% 
% fontsize = 14;
% figure;
% plot(d, [f_SqExp; f_SqExpNorm; f_Mat52; f_Mat32; f_Exp; f_ExpNorm; f_Sinc], 'linewidth', 1.5);
% xlabel('Distance d', 'fontsize', fontsize);
% ylabel('RBF f(d)', 'fontsize', fontsize); 
% title('Radial Basis Functions', 'fontsize', fontsize + 1);
% set(gca, 'fontsize', fontsize - 1);
% grid on; axis tight;
% h_lg = legend({'SqExp', 'SqExpNorm', 'Mat52',  'Mat32', 'Exp', 'ExpNorm',  'Sinc'}, 'location', 'best');
% set(h_lg, 'fontsize', fontsize -1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare Matlab versus Mex

% rng(340);
% d = rand(95, 65) * 10;
% ell = 1;
% f_SqExp_mat        = rbf_val('SqExp', d, ell, false);
% f_SqExp_mex        = rbf_val_mex('SqExp', d, ell, false);
% err = norm(f_SqExp_mat - f_SqExp_mex)

arguments
    mode (1,:) char  = 'SqExp';
    d (:,:) double {mustBeNonnegative} = 0;
    ell (1,1) double {mustBePositive} = 1;
    enable_disp (1,1) logical = false;
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if enable_disp && coder.target('MATLAB') 
    fontsize = 14;

    figure;
    plot(d, f, 'linewidth', 1.5);

    xlabel('Distance d', 'fontsize', fontsize);
    ylabel('RBF f(d)', 'fontsize', fontsize); 
    title(mode, 'fontsize', fontsize + 1, 'interpreter', 'none');
    set(gca, 'fontsize', fontsize - 1);
    grid on; axis tight;
end