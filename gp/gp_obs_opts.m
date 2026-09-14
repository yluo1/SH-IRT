function options = gp_obs_opts(options)
%Get default GP observations options

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input and output

%obs.omega:             [N_S x 1]    Angular frequency
%obs.theta:             [N_S x 1]    Co-latitude [0, pi]
%obs.phi:               [N_S x 1]    Azimuth [0, 2 * pi)
%obs.T60:               [N_S x 1]    Decay time in seconds
%obs.log_noise_std:     [N_S x 1]    Noise std. of log(T60), (std. as fraction of T60)

arguments
    options.omega            (:,1) double {mustBeNonnegative} = [];
    options.theta            (:,1) double = [];
    options.phi              (:,1) double = [];
    options.T60              (:,1) double {mustBeNonnegative} = [];
    options.log_noise_std    (:,1) double {mustBeNonnegative} = [];
end

assert(numel(options.omega) == numel(options.theta), 'options.omega, options.theta size mismatch');
assert(numel(options.omega) == numel(options.phi), 'options.omega, options.phi size mismatch');
assert(numel(options.omega) == numel(options.T60), 'options.omega, options.T60 size mismatch');
assert(numel(options.omega) == numel(options.log_noise_std), 'options.omega, options.log_noise_std size mismatch');