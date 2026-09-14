function options = gp_mu_opts(options)
%Get default GP prior mean options

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input and output

arguments
    options.mu_func    (1,:) char {mustBeMember(options.mu_func, {'power', 'LPF'}) } = 'power';

    options.mu_alpha   (1,1) double {mustBeNonnegative} = 1;
    options.mu_beta    (1,1) double {mustBeNonnegative} = 0.25;

    options.mu_fc      (1,1) double {mustBeNonnegative} = 8000;

    options.enable_disp (1,1) logical = false;
end