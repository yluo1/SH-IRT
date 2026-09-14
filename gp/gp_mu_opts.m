function options = gp_mu_opts(options)
%Get default GP prior mean options w.r.t. angular frequency

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input and output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    %Default covariance function
    options.mu_func    (1,:) char {mustBeMember(options.mu_func, {'power', 'LPF'}) } = 'power'; %Prior mean function name

    options.mu_alpha   (1,1) double {mustBeNonnegative} = 1;            % Normative function value
    options.mu_beta    (1,1) double {mustBeNonnegative} = 0.25;         % Decay exponent

    options.mu_fc      (1,1) double {mustBeNonnegative} = 8000;         % Cross-over frequency

    %Misc.
    options.enable_disp (1,1) logical = false;                          % Logical, if true, display plot in gp_mu.m
end