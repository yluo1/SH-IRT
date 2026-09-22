function [f, df_dell] = rbf_val(mode, d, ell, enable_disp)
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
%df_dell:       [N x M] Derivative of f w.r.t. ell

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
    mode (1,:) char {mustBeMember(mode, {'SqExp', 'SqExpNorm', 'Mat52', 'Mat32', 'Exp', 'ExpNorm', 'Sinc'})}  = 'SqExp';
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

    d_sq = d.^2;
    exp_term = exp(-d_sq ./ (2 * ell.^2));
    f = c * exp_term;   

    if nargout > 1
        dexp_dell = exp_term .* d_sq / (ell.^3);

        if contains(mode, 'SqExpNorm')
            exp_tmp = exp(-2 / (ell^2));
            dc_dell = ((ell^2 + 2) * exp_tmp  - ell^2) / (pi * (ell^5) * (1 - exp_tmp )^2);
            df_dell =  dc_dell * exp_term + c * dexp_dell;
        else
            df_dell = dexp_dell;
        end
    end

elseif strcmp(mode, 'Mat52')

    exp_term = exp(-sqrt(5) .* d ./ ell);
    f = (1 + sqrt(5) .* d ./ ell + 5/3 .* (d ./ ell).^2 ) .* exp_term;

    if nargout > 1
        df_dell = 5 * (d.^2) .* (ell + sqrt(5) * d) / (3 * ell^4) .* exp_term;
    end

elseif strcmp(mode, 'Mat32')

    exp_term = exp(-sqrt(3) * d / ell );
    f = (1 + sqrt(3) .* d ./ ell) .* exp_term; 
    if nargout > 1
        df_dell = 3 * (d.^2) / (ell^3) .* exp_term;
    end

elseif contains(mode, 'Exp')

    if contains(mode, 'ExpNorm')
        c =  1 / (2 * pi * (ell^2 - ell * exp(-2 / ell) * (ell + 2)));
    else
        c = 1;
    end
    
    exp_term = exp(-(d ./ ell));
    f = c * exp_term;

    if nargout > 1
        dexp_dell = exp_term .* d / (ell.^2);

        if contains(mode, 'ExpNorm')
            exp_tmp = exp(-2 / ell);
            dc_dell = (exp_tmp * (ell + 2 + 2 / ell) - ell) / (pi * (ell^2 - ell * (ell + 2) * exp_tmp )^2 );
            df_dell =  dc_dell * exp_term + c * dexp_dell;
        else
            df_dell = dexp_dell;
        end
    end
    
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